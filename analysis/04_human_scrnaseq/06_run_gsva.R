#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/06_run_gsva.R
# GSVA on single-cell data using pseudo-bulk approach
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/gsva_pseudobulk.rds
#   - analysis/04_human_scrnaseq/outputs/gsva_heatmap.pdf
#   - analysis/04_human_scrnaseq/figures/gsva_heatmap.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(GSVA)
  library(msigdbr)
  library(dplyr)
  library(data.table)
  library(pheatmap)
  library(tibble)
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
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== scRNA-seq GSVA (Pseudo-bulk) ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))
cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Create pseudo-bulk per patient
# -----------------------------------------------------------------------------
cat("\nStep 2: Creating pseudo-bulk profiles...\n")

# Aggregate counts per patient
DefaultAssay(seurat_obj) <- "RNA"

# Get raw counts
counts <- GetAssayData(seurat_obj, slot = "counts")

# Aggregate by patient (orig.ident)
patients <- unique(seurat_obj$orig.ident)
pseudobulk <- sapply(patients, function(pt) {
  cells <- colnames(seurat_obj)[seurat_obj$orig.ident == pt]
  rowSums(counts[, cells, drop = FALSE])
})

# Normalize (CPM + log)
pseudobulk_cpm <- sweep(pseudobulk, 2, colSums(pseudobulk), "/") * 1e6
pseudobulk_log <- log2(pseudobulk_cpm + 1)

cat("  Pseudo-bulk matrix:", nrow(pseudobulk_log), "genes x", ncol(pseudobulk_log), "samples\n")

# -----------------------------------------------------------------------------
# Step 3: Load gene sets (same as bulk analysis)
# -----------------------------------------------------------------------------
cat("\nStep 3: Loading gene sets...\n")

# Hallmark estrogen pathways
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
hallmark_estrogen <- hallmark_list[c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)]

# Reactome estrogen pathway
reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
reactome_list <- split(reactome_sets$gene_symbol, reactome_sets$gs_name)
reactome_estrogen <- reactome_list["REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION"]

# GO BP estrogen pathways
gobp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
gobp_list <- split(gobp_sets$gene_symbol, gobp_sets$gs_name)
gobp_estrogen <- gobp_list[c(
  "GOBP_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS"
)]

# Combine
estrogen_pathways <- c(hallmark_estrogen, reactome_estrogen, gobp_estrogen)
estrogen_pathways <- estrogen_pathways[!sapply(estrogen_pathways, is.null)]
cat("  Gene sets:", length(estrogen_pathways), "\n")

# -----------------------------------------------------------------------------
# Step 4: Run GSVA
# -----------------------------------------------------------------------------
cat("\nStep 4: Running GSVA...\n")

gsva_result <- gsva(
  gsvaParam(
    as.matrix(pseudobulk_log),
    estrogen_pathways,
    kcdf = "Gaussian",
    maxDiff = TRUE
  )
)

cat("  Result:", nrow(gsva_result), "pathways x", ncol(gsva_result), "samples\n")

# -----------------------------------------------------------------------------
# Step 5: Add age group annotation and visualize
# -----------------------------------------------------------------------------
cat("\nStep 5: Visualizing...\n")

# Get age group per patient
age_groups <- seurat_obj@meta.data %>%
  select(orig.ident, AgeGroup) %>%
  distinct()
rownames(age_groups) <- age_groups$orig.ident

# Annotation for heatmap
ann_col <- data.frame(
  AgeGroup = age_groups[colnames(gsva_result), "AgeGroup"],
  row.names = colnames(gsva_result)
)

ann_colors <- list(
  AgeGroup = c(Young = "#4DAF4A", MidAge = "#377EB8", Elderly = "#E41A1C")
)

pdf(file.path(output_dir, "gsva_heatmap.pdf"), width = 10, height = 6)
pheatmap(
  gsva_result,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  main = "GSVA Estrogen Pathway Scores (Pseudo-bulk)",
  scale = "row",
  show_colnames = TRUE
)
dev.off()

# Save PNG for validation pipeline (pheatmap needs explicit device)
png(file.path(figures_dir, "gsva_heatmap.png"), width = 10*300, height = 6*300, res = 300)
pheatmap(
  gsva_result,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  main = "GSVA Estrogen Pathway Scores (Pseudo-bulk)",
  scale = "row",
  show_colnames = TRUE
)
dev.off()

# -----------------------------------------------------------------------------
# Step 6: Save
# -----------------------------------------------------------------------------
cat("\nStep 6: Saving...\n")

saveRDS(gsva_result, file.path(output_dir, "gsva_pseudobulk.rds"))
saveRDS(pseudobulk_log, file.path(output_dir, "pseudobulk_log2cpm.rds"))

cat("\n=== scRNA GSVA complete ===\n")
