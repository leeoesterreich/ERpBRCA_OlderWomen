#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/15c_fig7d_elderly_minus_young.R
# Figure 7D - Elderly - Young HALLMARK GSVA z-score difference per cell type
# with row dendrogram + category column sidebar
# Reads pre-computed stats from 15_multicelltype_pathway.R
# Output: figures/fig7d_elderly_minus_young.{svg,pdf,png}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(showtext)
})

tryCatch(
  font_add("Arial", "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"),
  error = function(e) message("Arial not found; using default sans font")
)
showtext_auto()

# -----------------------------------------------------------------------------
# Load stats CSV (one row per celltype x pathway with `diff` column)
# -----------------------------------------------------------------------------
output_dir <- "analysis/04_human_scrnaseq/outputs"
figures_dir <- "analysis/04_human_scrnaseq/figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

stats_df <- fread(file.path(output_dir, "multicelltype_pathway_stats.csv"))

# Keep HALLMARK pathways only
stats_h <- stats_df[grepl("^HALLMARK_", pathway)]

# Pivot to cell type (cols) x pathway (rows)
mat <- stats_h |>
  dplyr::select(celltype, pathway, diff) |>
  tidyr::pivot_wider(names_from = celltype, values_from = diff) |>
  as.data.frame()
rownames(mat) <- mat$pathway
mat$pathway <- NULL
mat <- as.matrix(mat)

# Padj lookup for significance stars
sig_lookup <- setNames(stats_h$padj, paste0(stats_h$celltype, "::", stats_h$pathway))
sig_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
for (i in seq_len(nrow(mat))) {
  pw <- rownames(mat)[i]
  for (j in seq_len(ncol(mat))) {
    ct <- colnames(mat)[j]
    key <- paste0(ct, "::", pw)
    padj <- sig_lookup[key]
    if (!is.na(padj) && padj < 0.05) sig_mat[i, j] <- "*"
    if (!is.na(padj) && padj < 0.01) sig_mat[i, j] <- "**"
  }
}

# -----------------------------------------------------------------------------
# Category map for columns (same as 7B)
# -----------------------------------------------------------------------------
get_category <- function(ct) {
  dplyr::case_when(
    grepl("Tcells|NK|Bcells|Plasma|Cycling", ct) & !grepl("Myeloid", ct) ~ "Lymphocyte",
    grepl("Macro|Mono|DC|Myeloid", ct) ~ "Myeloid",
    grepl("Epithelial|Cancer", ct) ~ "Epithelial",
    grepl("CAF|PVL|Endo", ct) ~ "Stromal",
    TRUE ~ "Other"
  )
}
category_order <- c("Lymphocyte", "Myeloid", "Epithelial", "Stromal", "Other")
col_categories <- sapply(colnames(mat), get_category)
col_order <- order(match(col_categories, category_order), colnames(mat))
mat <- mat[, col_order, drop = FALSE]
sig_mat <- sig_mat[, col_order, drop = FALSE]
col_categories <- col_categories[col_order]

# Display rownames: strip HALLMARK_ prefix
rownames(mat) <- gsub("^HALLMARK_", "", rownames(mat))
rownames(mat) <- gsub("_", " ", rownames(mat))

ann_colors <- list(
  Category = c(Lymphocyte = "#009E73", Myeloid = "#D55E00",
               Epithelial = "#0072B2", Stromal = "#CC79A7", Other = "#000000")
)

top_anno <- HeatmapAnnotation(
  Category = col_categories,
  col = ann_colors,
  show_annotation_name = FALSE,
  height = unit(3, "mm"),
  border = FALSE,
  which = "column"
)

# Diff range is roughly ±0.8; use ±0.6 as in PDF
diff_range <- 0.6
diff_col_fun <- colorRamp2(
  c(-diff_range, -diff_range/2, 0, diff_range/2, diff_range),
  c("#3E5CA8", "#8BC5E4", "#EDF3F1", "#F6BD6A", "#E75321")
)

ht <- Heatmap(
  mat,
  name = "Elderly - Young",
  col = diff_col_fun,
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  row_dend_width = unit(12, "mm"),
  column_split = factor(col_categories, levels = category_order),
  column_gap = unit(2, "mm"),
  column_title_gp = gpar(fontsize = 11, fontfamily = "Arial"),
  top_annotation = top_anno,
  show_row_names = TRUE, row_names_side = "right",
  row_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  show_column_names = TRUE, column_names_rot = 45,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  column_title = "Pathway Activity: Elderly - Young (GSVA z-score difference)",
  border = FALSE,
  rect_gp = gpar(col = "#3B3B3B", lwd = 0.4),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] != "") {
      grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 8, col = "black", fontfamily = "Arial"))
    }
  },
  width = unit(100, "mm"),
  height = unit(140, "mm"),
  heatmap_legend_param = list(
    direction = "vertical",
    title = "Elderly - Young\n(GSVA diff)",
    at = c(-diff_range, 0, diff_range),
    legend_height = unit(35, "mm"),
    title_gp = gpar(fontsize = 10, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 9, fontfamily = "Arial")
  )
)

svg_path <- file.path(figures_dir, "fig7d_elderly_minus_young.svg")
pdf_path <- file.path(figures_dir, "fig7d_elderly_minus_young.pdf")
png_path <- file.path(figures_dir, "fig7d_elderly_minus_young.png")

svg(svg_path, width = 8.5, height = 8.0)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
invisible(dev.off())
cat("Saved:", svg_path, "\n")

pdf(pdf_path, width = 8.5, height = 8.0)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
invisible(dev.off())
cat("Saved:", pdf_path, "\n")

png(png_path, width = 8.5, height = 8.0, units = "in", res = 300)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
invisible(dev.off())
cat("Saved:", png_path, "\n")
