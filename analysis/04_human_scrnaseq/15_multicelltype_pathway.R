#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/15_multicelltype_pathway.R
# HALLMARK/BIOCARTA pathway enrichment per cell type
# Generates Figure 7 Panels B/C - dual pathway activity heatmaps
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/macrophage_seurat.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/multicelltype_pathway_scores.csv
#   - analysis/04_human_scrnaseq/figures/fig7bc_pathway_heatmaps.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(msigdbr)
  library(GSVA)
  library(pheatmap)
  library(tibble)
  library(gridExtra)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "outputs")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Multi-Cell-Type Pathway Enrichment (Figure 7C) ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_file <- file.path(output_dir, "macrophage_seurat.rds")
seurat_obj <- readRDS(seurat_file)

# Filter to Elderly and Young
seurat_ey <- subset(seurat_obj, subset = AgeGroup %in% c("Elderly", "Young"))
cat("  Cells:", ncol(seurat_ey), "\n")

# -----------------------------------------------------------------------------
# Step 2: Load HALLMARK and BIOCARTA gene sets
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading gene sets...\n")

# HALLMARK - use ALL pathways (50 total) for better z-score range
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
cat("  HALLMARK pathways loaded:", length(hallmark_list), "\n")

# BIOCARTA - use ALL pathways
biocarta_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
biocarta_list <- split(biocarta_sets$gene_symbol, biocarta_sets$gs_name)
cat("  BIOCARTA pathways loaded:", length(biocarta_list), "\n")

all_pathways <- c(hallmark_list, biocarta_list)
cat("  Total pathways:", length(all_pathways), "\n")

# -----------------------------------------------------------------------------
# Step 3: Create pseudo-bulk per cell type per age group
# -----------------------------------------------------------------------------
cat("\nStep 3: Creating pseudo-bulk profiles...\n")

DefaultAssay(seurat_ey) <- "RNA"
counts <- GetAssayData(seurat_ey, slot = "counts")

# Create grouping variable
seurat_ey$ct_age <- paste0(seurat_ey$CellTypeAnnotSH, "_", seurat_ey$AgeGroup)
groups <- unique(seurat_ey$ct_age)

# Aggregate
pseudobulk <- sapply(groups, function(grp) {
  cells <- colnames(seurat_ey)[seurat_ey$ct_age == grp]
  if (length(cells) > 10) {
    rowSums(counts[, cells, drop = FALSE])
  } else {
    rep(NA, nrow(counts))
  }
})

# Remove NA columns
pseudobulk <- pseudobulk[, !apply(pseudobulk, 2, function(x) all(is.na(x)))]

# Normalize
pseudobulk_cpm <- sweep(pseudobulk, 2, colSums(pseudobulk, na.rm = TRUE), "/") * 1e6
pseudobulk_log <- log2(pseudobulk_cpm + 1)

cat("  Pseudo-bulk matrix:", nrow(pseudobulk_log), "genes x", ncol(pseudobulk_log), "groups\n")

# -----------------------------------------------------------------------------
# Step 4: Run GSVA
# -----------------------------------------------------------------------------
cat("\nStep 4: Running GSVA...\n")

gsva_result <- gsva(
  gsvaParam(
    as.matrix(pseudobulk_log),
    all_pathways,
    kcdf = "Gaussian",
    maxDiff = TRUE
  )
)

cat("  GSVA result:", nrow(gsva_result), "pathways x", ncol(gsva_result), "groups\n")

# Per-pathway z-score normalization (row-wise) to expand range
cat("  Applying per-pathway z-score normalization...\n")
gsva_result <- t(scale(t(gsva_result)))
cat("  Z-score range:", round(min(gsva_result, na.rm = TRUE), 2), "to",
    round(max(gsva_result, na.rm = TRUE), 2), "\n")

# Save normalized scores
gsva_df <- as.data.frame(gsva_result) %>%
  tibble::rownames_to_column("pathway")
fwrite(gsva_df, file.path(output_dir, "multicelltype_pathway_scores.csv"))
cat("  Saved multicelltype_pathway_scores.csv\n")

# -----------------------------------------------------------------------------
# Step 5: Generate Figure 7 B/C - Dual Heatmaps (HALLMARK + BIOCARTA)
# -----------------------------------------------------------------------------
cat("\nStep 5: Generating Figure 7 B/C dual heatmaps...\n")

# Transpose: rows = cell_type_age, columns = pathways
gsva_t <- t(gsva_result)

# Create row annotations for cell type categories
get_category <- function(ct_age) {
  ct <- gsub("_(Elderly|Young)$", "", ct_age)
  case_when(
    grepl("Tcells|NK|Bcells|Plasma", ct) ~ "Lymphocyte",
    grepl("Macro|Mono|DC|Myeloid", ct) ~ "Myeloid",
    grepl("Epithelial|Cancer|Luminal|Basal", ct) ~ "Epithelial",
    grepl("CAF|PVL|Endo|Fibro", ct) ~ "Stromal",
    TRUE ~ "Other"
  )
}

row_categories <- data.frame(
  Category = sapply(rownames(gsva_t), get_category),
  row.names = rownames(gsva_t)
)

# Define row order: group by category, then by cell type (Younger before Older)
category_order <- c("Lymphocyte", "Myeloid", "Epithelial", "Stromal", "Other")
row_order <- rownames(gsva_t)[order(
  match(row_categories$Category, category_order),
  gsub("_(Elderly|Young)$", "", rownames(gsva_t)),
  grepl("_Elderly$", rownames(gsva_t))
)]
gsva_t <- gsva_t[row_order, ]
row_categories <- row_categories[row_order, , drop = FALSE]

# Calculate gaps for category separation
category_counts <- table(row_categories$Category)[category_order]
category_counts <- category_counts[!is.na(category_counts) & category_counts > 0]
gaps_row <- cumsum(category_counts)[-length(category_counts)]

# Split into HALLMARK and BIOCARTA
hallmark_cols <- grep("^HALLMARK_", colnames(gsva_t), value = TRUE)
biocarta_cols <- grep("^BIOCARTA_", colnames(gsva_t), value = TRUE)

gsva_hallmark <- gsva_t[, hallmark_cols, drop = FALSE]
gsva_biocarta <- gsva_t[, biocarta_cols, drop = FALSE]

# Clean column names for display
colnames(gsva_hallmark) <- gsub("^HALLMARK_", "", colnames(gsva_hallmark))
colnames(gsva_hallmark) <- gsub("_", " ", colnames(gsva_hallmark))
colnames(gsva_biocarta) <- gsub("^BIOCARTA_", "", colnames(gsva_biocarta))
colnames(gsva_biocarta) <- gsub("_", " ", colnames(gsva_biocarta))

# Annotation colors
ann_colors <- list(
  Category = c(
    Lymphocyte = "#4DAF4A",
    Myeloid = "#E41A1C",
    Epithelial = "#377EB8",
    Stromal = "#984EA3",
    Other = "grey70"
  )
)

# Color scale (blue-white-red, symmetric around 0)
max_val <- max(abs(gsva_t), na.rm = TRUE)
color_breaks <- seq(-max_val, max_val, length.out = 101)
color_palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Create HALLMARK heatmap
p_hallmark <- pheatmap(
  gsva_hallmark,
  annotation_row = row_categories,
  annotation_colors = ann_colors,
  main = "HALLMARK",
  color = color_palette,
  breaks = color_breaks,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  gaps_row = gaps_row,
  fontsize_row = 8,
  fontsize_col = 9,
  show_rownames = TRUE,
  angle_col = 45,
  silent = TRUE
)

# Create BIOCARTA heatmap
p_biocarta <- pheatmap(
  gsva_biocarta,
  annotation_row = row_categories,
  annotation_colors = ann_colors,
  main = "BIOCARTA",
  color = color_palette,
  breaks = color_breaks,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  gaps_row = gaps_row,
  fontsize_row = 8,
  fontsize_col = 9,
  show_rownames = TRUE,
  angle_col = 45,
  silent = TRUE
)

# Combine side-by-side
png(file.path(figures_dir, "fig7bc_pathway_heatmaps.png"),
    width = 16*300, height = 12*300, res = 300)
grid.arrange(p_hallmark$gtable, p_biocarta$gtable, ncol = 2)
dev.off()
cat("  Saved fig7bc_pathway_heatmaps.png\n")

cat("\n=== Multi-cell-type pathway analysis complete ===\n")
