#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/02_harmony_integrate.R
# Harmony batch correction across patients
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_young_midage_elderly.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_harmony_integrated.rds
#   - analysis/04_human_scrnaseq/outputs/umap_by_patient.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
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

cat("=== Harmony Integration ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_young_midage_elderly.rds"))
cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Run Harmony
# -----------------------------------------------------------------------------
cat("\nStep 2: Running Harmony integration...\n")

seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "orig.ident",
  reduction = "pca",
  dims.use = 1:30,
  verbose = FALSE
)

# Run UMAP on Harmony reduction
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.5)

cat("  Clusters found:", length(unique(seurat_obj$seurat_clusters)), "\n")

# -----------------------------------------------------------------------------
# Step 3: Visualize
# -----------------------------------------------------------------------------
cat("\nStep 3: Generating visualizations...\n")

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
  ggtitle("By Patient")

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "AgeGroup") +
  ggtitle("By Age Group")

p3 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
  ggtitle("Clusters")

pdf(file.path(output_dir, "umap_harmony.pdf"), width = 15, height = 5)
print(p1 + p2 + p3)
dev.off()

# -----------------------------------------------------------------------------
# Step 4: Save
# -----------------------------------------------------------------------------
cat("\nStep 4: Saving...\n")

saveRDS(seurat_obj, file.path(output_dir, "seurat_harmony_integrated.rds"))

cat("\n=== Harmony complete ===\n")
