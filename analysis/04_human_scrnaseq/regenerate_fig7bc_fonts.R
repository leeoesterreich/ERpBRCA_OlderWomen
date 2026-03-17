#!/usr/bin/env Rscript
# Regenerate Fig 7B/C heatmaps from existing GSVA scores
# Only changes: Arial font + increased font sizes
# Reads from: outputs/multicelltype_pathway_scores.csv, outputs/multicelltype_pathway_stats.csv
# Outputs: figures/fig7bc_pathway_heatmaps.png/.svg

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("--file=", "", file_arg))))
  return(getwd())
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "outputs")
figures_dir <- file.path(script_dir, "figures")

cat("=== Regenerating Fig 7B/C with Arial font + larger sizes ===\n")

# Load pre-computed GSVA z-scores
scores <- fread(file.path(output_dir, "multicelltype_pathway_scores.csv"))
gsva_z <- as.matrix(scores[, -1, with = FALSE])
rownames(gsva_z) <- scores$celltype_age

# Load stats for significance stars
stats_df <- fread(file.path(output_dir, "multicelltype_pathway_stats.csv"))

cat("  Loaded", nrow(gsva_z), "rows x", ncol(gsva_z), "pathways\n")

# Row categories
get_category <- function(ct_age) {
  ct <- gsub("_(Elderly|Young|Younger|Older)$", "", ct_age)
  dplyr::case_when(
    grepl("Tcells|NK|Bcells|Plasma|Cycling", ct) & !grepl("Myeloid", ct) ~ "Lymphocyte",
    grepl("Macro|Mono|DC|Myeloid", ct) ~ "Myeloid",
    grepl("Epithelial|Cancer", ct) ~ "Epithelial",
    grepl("CAF|PVL|Endo", ct) ~ "Stromal",
    TRUE ~ "Other"
  )
}

category_order <- c("Lymphocyte", "Myeloid", "Epithelial", "Stromal")
row_categories <- data.frame(
  Category = sapply(rownames(gsva_z), get_category),
  row.names = rownames(gsva_z)
)

row_order <- rownames(gsva_z)[order(
  match(row_categories$Category, category_order),
  gsub("_(Younger|Older)$", "", rownames(gsva_z)),
  grepl("_Older$", rownames(gsva_z))
)]
gsva_z <- gsva_z[row_order, ]
row_categories <- row_categories[row_order, , drop = FALSE]

# Select pathways
hallmark_cols <- grep("^(HALLMARK_|LI_ESTROGENE_)", colnames(gsva_z), value = TRUE)
biocarta_cols <- grep("^BIOCARTA_", colnames(gsva_z), value = TRUE)

hallmark_curated <- c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS",
  "LI_ESTROGENE_EARLY_E2_RESPONSE_UP", "LI_ESTROGENE_LATE_E2_RESPONSE_UP"
)
biocarta_curated <- c(
  "BIOCARTA_INFLAM_PATHWAY", "BIOCARTA_IL6_PATHWAY", "BIOCARTA_IL2_PATHWAY",
  "BIOCARTA_NFKB_PATHWAY", "BIOCARTA_TNFR1_PATHWAY", "BIOCARTA_DEATH_PATHWAY",
  "BIOCARTA_FAS_PATHWAY", "BIOCARTA_CASPASE_PATHWAY", "BIOCARTA_P53_PATHWAY",
  "BIOCARTA_CELLCYCLE_PATHWAY", "BIOCARTA_G1_PATHWAY", "BIOCARTA_G2_PATHWAY"
)
hallmark_keep <- intersect(hallmark_curated, hallmark_cols)
biocarta_keep <- intersect(biocarta_curated, biocarta_cols)
cat("  Using", length(hallmark_keep), "HALLMARK and", length(biocarta_keep), "BIOCARTA pathways\n")

gsva_hallmark <- gsva_z[, hallmark_keep, drop = FALSE]
gsva_biocarta <- gsva_z[, biocarta_keep, drop = FALSE]

# Significance stars
sig_lookup <- setNames(stats_df$padj, paste0(stats_df$celltype, "::", stats_df$pathway))

make_sig_matrix <- function(mat, pw_names_original) {
  sig_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
  for (i in 1:nrow(mat)) {
    rn <- rownames(mat)[i]
    ct <- gsub("_(Younger|Older)$", "", rn)
    for (j in 1:ncol(mat)) {
      key <- paste0(ct, "::", pw_names_original[j])
      padj <- sig_lookup[key]
      if (!is.na(padj) && padj < 0.05) sig_mat[i, j] <- "*"
      if (!is.na(padj) && padj < 0.01) sig_mat[i, j] <- "**"
    }
  }
  sig_mat
}

hallmark_sig <- make_sig_matrix(gsva_hallmark, hallmark_keep)
biocarta_sig <- make_sig_matrix(gsva_biocarta, biocarta_keep)

# Clean column names
colnames(gsva_hallmark) <- gsub("^(HALLMARK_|LI_ESTROGENE_)", "", colnames(gsva_hallmark))
colnames(gsva_hallmark) <- gsub("_", " ", colnames(gsva_hallmark))
colnames(gsva_biocarta) <- gsub("^BIOCARTA_", "", colnames(gsva_biocarta))
colnames(gsva_biocarta) <- gsub("_", " ", colnames(gsva_biocarta))

# Heatmap setup
ann_colors <- list(
  Category = c(Lymphocyte = "#4DAF4A", Myeloid = "#E41A1C", Epithelial = "#377EB8", Stromal = "#984EA3")
)

hallmark_col_fun <- colorRamp2(c(-3, -1.5, 0, 1.5, 3), c("#3E5CA8", "#8BC5E4", "#EDF3F1", "#F6BD6A", "#E75321"))
biocarta_col_fun <- colorRamp2(c(-4, -2, 0, 2, 4), c("#3E5CA8", "#8BC5E4", "#EDF3F1", "#F6BD6A", "#E75321"))

row_split <- factor(row_categories$Category, levels = category_order)

left_anno <- rowAnnotation(
  Category = row_categories$Category,
  col = ann_colors,
  show_annotation_name = FALSE,
  width = unit(3, "mm"),
  border = FALSE
)

make_heatmap <- function(mat, sig_mat, title_text, col_fun, range_vals) {
  Heatmap(
    mat, name = "z-score", col = col_fun,
    cluster_rows = FALSE, cluster_columns = TRUE,
    show_row_dend = FALSE,
    row_split = row_split, row_gap = unit(2, "mm"),
    row_title_side = "left", row_title_rot = 0,
    row_title_gp = gpar(fontsize = 12, fontface = "plain", fontfamily = "Arial"),
    left_annotation = left_anno,
    show_row_names = TRUE, row_names_side = "right",
    row_names_gp = gpar(fontsize = 11, fontfamily = "Arial"),
    show_column_names = TRUE, column_names_rot = 60,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
    column_title = title_text,
    column_title_gp = gpar(fontsize = 18, fontface = "plain", fontfamily = "Arial"),
    column_dend_height = unit(10, "mm"),
    border = FALSE,
    rect_gp = gpar(col = "#3B3B3B", lwd = 0.5),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (sig_mat[i, j] != "") {
        grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 10, col = "black", fontfamily = "Arial"))
      }
    },
    width = unit(55, "mm"),
    height = unit(110, "mm"),
    heatmap_legend_param = list(
      direction = "horizontal",
      title = "Pathway activity\n(z-score)",
      at = c(range_vals[1], 0, range_vals[2]),
      legend_width = unit(30, "mm"),
      title_gp = gpar(fontsize = 12, fontfamily = "Arial"),
      labels_gp = gpar(fontsize = 11, fontfamily = "Arial")
    )
  )
}

ht_hallmark <- make_heatmap(gsva_hallmark, hallmark_sig, "Hallmark", hallmark_col_fun, c(-3, 3))
ht_biocarta <- make_heatmap(gsva_biocarta, biocarta_sig, "BIOCARTA", biocarta_col_fun, c(-4, 4))

png_file <- file.path(figures_dir, "fig7bc_pathway_heatmaps.png")
png(png_file, width = 3000, height = 2400, res = 300, bg = "white")
draw(
  ht_hallmark + ht_biocarta,
  heatmap_legend_side = "bottom",
  merge_legends = TRUE,
  padding = unit(c(10, 15, 22, 10), "mm")
)
dev.off()
cat("  Saved", basename(png_file), "\n")

svg_file <- file.path(figures_dir, "fig7bc_pathway_heatmaps.svg")
svg(svg_file, width = 10, height = 8)
draw(
  ht_hallmark + ht_biocarta,
  heatmap_legend_side = "bottom",
  merge_legends = TRUE,
  padding = unit(c(10, 15, 22, 10), "mm")
)
dev.off()
cat("  Saved", basename(svg_file), "\n")

cat("\n=== Done ===\n")
