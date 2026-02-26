#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/07_run_progeny.R
# PROGENy pathway activity on single-cell data
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/progeny_activity.rds
#   - analysis/04_human_scrnaseq/outputs/progeny_heatmap.pdf
#   - analysis/04_human_scrnaseq/figures/progeny_heatmap.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(progeny)
  library(dplyr)
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
output_dir <- file.path(script_dir, "outputs")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== scRNA-seq PROGENy ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))

# Use normalized data
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Create pseudo-bulk per patient
# -----------------------------------------------------------------------------
cat("\nStep 2: Creating pseudo-bulk profiles...\n")

# Use the log-normalized data
norm_data <- GetAssayData(seurat_obj, slot = "data")

# Aggregate by patient (mean expression)
patients <- unique(seurat_obj$orig.ident)
pseudobulk <- sapply(patients, function(pt) {
  cells <- colnames(seurat_obj)[seurat_obj$orig.ident == pt]
  rowMeans(norm_data[, cells, drop = FALSE])
})

cat("  Pseudo-bulk matrix:", nrow(pseudobulk), "genes x", ncol(pseudobulk), "samples\n")

# -----------------------------------------------------------------------------
# Step 3: Run PROGENy
# -----------------------------------------------------------------------------
cat("\nStep 3: Running PROGENy...\n")

pathway_activity <- progeny(
  as.matrix(pseudobulk),
  scale = FALSE,
  organism = "Human",
  top = 100,
  perm = 1000
)

cat("  Result:", nrow(pathway_activity), "samples x", ncol(pathway_activity), "pathways\n")

# -----------------------------------------------------------------------------
# Step 4: Visualize
# -----------------------------------------------------------------------------
cat("\nStep 4: Visualizing...\n")

# Get age group per patient
age_groups <- seurat_obj@meta.data %>%
  select(orig.ident, AgeGroup) %>%
  distinct()
rownames(age_groups) <- age_groups$orig.ident

# Annotation
ann_row <- data.frame(
  AgeGroup = age_groups[rownames(pathway_activity), "AgeGroup"],
  row.names = rownames(pathway_activity)
)

ann_colors <- list(
  AgeGroup = c(Young = "#4DAF4A", MidAge = "#377EB8", Elderly = "#E41A1C")
)

pdf(file.path(output_dir, "progeny_heatmap.pdf"), width = 10, height = 8)
pheatmap(
  pathway_activity,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  main = "PROGENy Pathway Activity (Pseudo-bulk)",
  scale = "column",
  cluster_cols = FALSE
)
dev.off()

# Save PNG for validation pipeline (pheatmap needs explicit device)
png(file.path(figures_dir, "progeny_heatmap.png"), width = 10*300, height = 8*300, res = 300)
pheatmap(
  pathway_activity,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  main = "PROGENy Pathway Activity (Pseudo-bulk)",
  scale = "column",
  cluster_cols = FALSE
)
dev.off()

# Focus on Estrogen pathway
if ("Estrogen" %in% colnames(pathway_activity)) {
  estrogen_activity <- data.frame(
    Patient = rownames(pathway_activity),
    Estrogen = pathway_activity[, "Estrogen"],
    AgeGroup = ann_row$AgeGroup
  )
  write.csv(estrogen_activity, file.path(output_dir, "progeny_estrogen_activity.csv"),
            row.names = FALSE)
}

# -----------------------------------------------------------------------------
# Step 5: Save
# -----------------------------------------------------------------------------
cat("\nStep 5: Saving...\n")

saveRDS(pathway_activity, file.path(output_dir, "progeny_activity.rds"))

cat("\n=== scRNA PROGENy complete ===\n")
