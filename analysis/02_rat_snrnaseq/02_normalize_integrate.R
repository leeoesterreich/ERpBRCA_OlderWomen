#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/02_normalize_integrate.R
# SCTransform normalization and Harmony integration for rat snRNA-seq
#
# Inputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_qc_filtered.rds
#   - data/rat_snrnaseq/metadata/Rat_scRNAseq_AgeGroup.txt
#
# Outputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_integrated.rds
#   - analysis/02_rat_snrnaseq/outputs/integration_plots.pdf

# Set random seed for reproducibility
set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(ggplot2)
  library(dplyr)
})

# Helper function to check file existence
check_file_exists <- function(filepath, description = "file") {
  if (!file.exists(filepath)) {
    stop(sprintf("ERROR: %s not found: %s", description, filepath))
  }
  cat(sprintf("  Found: %s\n", basename(filepath)))
}

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
metadata_dir <- file.path(project_root, "data/rat_snrnaseq/metadata")

cat("=== Rat snRNA-seq Normalization and Integration ===\n")
cat("Random seed: 12345\n")
cat("Project root:", project_root, "\n")
cat("Output directory:", output_dir, "\n\n")

# Input files
input_rds <- file.path(output_dir, "seurat_qc_filtered.rds")
age_group_file <- file.path(metadata_dir, "Rat_scRNAseq_AgeGroup.txt")

# Verify input files exist
cat("Checking input files...\n")
check_file_exists(input_rds, "QC-filtered Seurat object")
check_file_exists(age_group_file, "Age group annotation file")
cat("\n")

# -----------------------------------------------------------------------------
# Step 1: Load Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading QC-filtered Seurat object...\n")
seurat_obj <- readRDS(input_rds)
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Samples:", length(unique(seurat_obj$orig.ident)), "\n")

# Load age group metadata
age_groups <- read.delim(age_group_file, stringsAsFactors = FALSE)
cat("  Age group annotations loaded:", nrow(age_groups), "samples\n")
print(age_groups)

# -----------------------------------------------------------------------------
# Step 2: Add Age Group Metadata
# -----------------------------------------------------------------------------
cat("\nStep 2: Adding age group metadata...\n")

# Map CaseID to AgeGroup
seurat_obj$AgeGroup <- age_groups$AgeGroup[match(seurat_obj$orig.ident, age_groups$CaseID)]

# Verify mapping
if (any(is.na(seurat_obj$AgeGroup))) {
  warning("Some cells have missing AgeGroup assignments")
  cat("  Samples with missing mapping:\n")
  print(unique(seurat_obj$orig.ident[is.na(seurat_obj$AgeGroup)]))
}

cat("  Age group distribution:\n")
print(table(seurat_obj$AgeGroup))

# -----------------------------------------------------------------------------
# Step 3: SCTransform Normalization on Merged Object
# -----------------------------------------------------------------------------
cat("\nStep 3: Re-running SCTransform on merged object...\n")

# Re-run SCTransform on the merged object
seurat_obj <- SCTransform(
  seurat_obj,
  method = "glmGamPoi",
  vars.to.regress = "percent.mt",
  verbose = TRUE,
  seed.use = 12345
)

cat("  SCTransform complete\n")

# -----------------------------------------------------------------------------
# Step 4: PCA
# -----------------------------------------------------------------------------
cat("\nStep 4: Running PCA...\n")
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE, seed.use = 12345)

# Determine optimal number of PCs using elbow plot data
elbow_data <- ElbowPlot(seurat_obj, ndims = 50, reduction = "pca")
cat("  PCA complete, 50 PCs computed\n")

# -----------------------------------------------------------------------------
# Step 5: Harmony Integration by Sample
# -----------------------------------------------------------------------------
cat("\nStep 5: Running Harmony integration...\n")

seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "orig.ident",
  reduction = "pca",
  assay.use = "SCT",
  plot_convergence = FALSE,
  seed = 12345
)

cat("  Harmony integration complete\n")

# -----------------------------------------------------------------------------
# Step 6: Generate Integration QC Plots
# -----------------------------------------------------------------------------
cat("\nStep 6: Generating integration QC plots...\n")

# Run UMAP on harmony embeddings for visualization
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  seed.use = 12345
)

pdf(file.path(output_dir, "integration_plots.pdf"), width = 14, height = 10)

# Elbow plot
print(ElbowPlot(seurat_obj, ndims = 50, reduction = "pca") +
  ggtitle("PCA Elbow Plot") +
  theme_bw())

# UMAP by sample (pre-integration view)
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP colored by Sample") +
  theme_bw()
print(p1)

# UMAP by age group
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "AgeGroup") +
  ggtitle("UMAP colored by Age Group") +
  theme_bw()
print(p2)

# Split UMAP by sample
p3 <- DimPlot(seurat_obj, reduction = "umap", split.by = "orig.ident", ncol = 3) +
  ggtitle("UMAP split by Sample") +
  theme_bw()
print(p3)

# Split UMAP by age group
p4 <- DimPlot(seurat_obj, reduction = "umap", split.by = "AgeGroup") +
  ggtitle("UMAP split by Age Group") +
  theme_bw()
print(p4)

# Feature plots for QC metrics
p5 <- FeaturePlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3) &
  theme_bw()
print(p5)

dev.off()

cat("  Integration plots saved to:", file.path(output_dir, "integration_plots.pdf"), "\n")

# -----------------------------------------------------------------------------
# Step 7: Save Integrated Object
# -----------------------------------------------------------------------------
cat("\nStep 7: Saving integrated object...\n")
saveRDS(seurat_obj, file.path(output_dir, "seurat_integrated.rds"))

# Memory cleanup
gc()

cat("\n=== Normalization and integration complete ===\n")
cat("Output:", file.path(output_dir, "seurat_integrated.rds"), "\n")
