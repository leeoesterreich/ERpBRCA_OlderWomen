#!/usr/bin/env Rscript
# Fast re-render of Fig 2E using pre-computed CSVs from 10_multicelltype_pathway.R
# Lets us iterate on layout without re-running the full GSVA pipeline.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(showtext)
})

tryCatch(font_add("Arial", "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"),
         error = function(e) message("Arial not found; using default sans font"))
showtext_auto()

output_dir <- "analysis/02_rat_snrnaseq/outputs"
figures_dir <- "analysis/02_rat_snrnaseq/figures"

scores <- fread(file.path(output_dir, "rat_multicelltype_pathway_scores.csv"))
stats_df <- fread(file.path(output_dir, "rat_multicelltype_pathway_stats.csv"))

gsva_z <- as.matrix(scores[, -1])
rownames(gsva_z) <- scores$celltype_age
mat <- t(gsva_z)

get_category <- function(ct) {
  dplyr::case_when(
    grepl("CD4Tcell|CD8Tcell|NaiveTcell|NKcell|Treg", ct) ~ "Lymphocyte",
    grepl("Macrophage|Monocyte|DendriticCell", ct) ~ "Myeloid",
    grepl("CancerEpithelial|Myoepithelial", ct) ~ "Epithelial",
    grepl("Endothelial|Fibroblast", ct) ~ "Stromal",
    TRUE ~ "Other"
  )
}

col_celltypes <- sub("_(Young|Elderly)$", "", colnames(mat))
col_categories <- sapply(col_celltypes, get_category)
category_order <- c("Lymphocyte", "Myeloid", "Epithelial", "Stromal", "Other")
col_order <- order(match(col_categories, category_order), col_celltypes,
                   grepl("_Elderly$", colnames(mat)))
mat <- mat[, col_order, drop = FALSE]
col_categories <- col_categories[col_order]

display_cols <- gsub("_Young$", "_Younger", colnames(mat))
display_cols <- gsub("_Elderly$", "_Older", display_cols)
colnames(mat) <- display_cols
display_rows <- gsub("^HALLMARK_", "", rownames(mat))
display_rows <- gsub("_", " ", display_rows)
rownames(mat) <- display_rows

sig_lookup <- setNames(stats_df$padj, paste0(stats_df$celltype, "::", stats_df$pathway))
sig_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
for (j in seq_len(ncol(mat))) {
  ct_age <- colnames(mat)[j]
  ct <- gsub("_(Younger|Older)$", "", ct_age)
  for (i in seq_len(nrow(mat))) {
    pw <- paste0("HALLMARK_", gsub(" ", "_", rownames(mat)[i]))
    padj <- sig_lookup[paste0(ct, "::", pw)]
    if (!is.na(padj) && padj < 0.05) sig_mat[i, j] <- "*"
    if (!is.na(padj) && padj < 0.01) sig_mat[i, j] <- "**"
  }
}

ann_colors <- list(Category = c(
  Lymphocyte = "#009E73", Myeloid = "#D55E00",
  Epithelial = "#0072B2", Stromal = "#CC79A7", Other = "#000000"))
col_split <- factor(col_categories, levels = category_order)
col_fun <- colorRamp2(c(-3, -1.5, 0, 1.5, 3),
                      c("#3E5CA8", "#8BC5E4", "#EDF3F1", "#F6BD6A", "#E75321"))

bottom_anno <- HeatmapAnnotation(
  Category = col_categories, col = ann_colors,
  show_annotation_name = FALSE, height = unit(3, "mm"),
  border = FALSE, which = "column"
)

ht <- Heatmap(
  mat, name = "Pathway activity", col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_split, column_gap = unit(2, "mm"),
  column_title_gp = gpar(fontsize = 11, fontfamily = "Arial"),
  bottom_annotation = bottom_anno,
  show_row_names = TRUE, row_names_side = "right",
  row_names_gp = gpar(fontsize = 10, fontfamily = "Arial"),
  show_column_names = TRUE, column_names_rot = 45,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  border = FALSE, rect_gp = gpar(col = "#3B3B3B", lwd = 0.4),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] != "") {
      grid.text(sig_mat[i, j], x, y,
                gp = gpar(fontsize = 8, col = "black", fontfamily = "Arial"))
    }
  },
  width = unit(140, "mm"),
  height = unit(190, "mm"),
  heatmap_legend_param = list(
    direction = "horizontal",
    title = "Pathway activity\n(z-score)",
    at = c(-3, 0, 3), legend_width = unit(30, "mm"),
    title_gp = gpar(fontsize = 10, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 9, fontfamily = "Arial")
  )
)

for (ext in c("svg", "pdf", "png")) {
  path <- file.path(figures_dir, paste0("fig2e_rat_pathway_heatmap.", ext))
  if (ext == "png") {
    png(path, width = 10, height = 10.5, units = "in", res = 300)
  } else if (ext == "pdf") {
    pdf(path, width = 10, height = 10.5)
  } else {
    svg(path, width = 10, height = 10.5)
  }
  draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
  invisible(dev.off())
  cat("Saved:", path, "\n")
}
