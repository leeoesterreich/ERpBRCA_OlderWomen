#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/00b_preprocess_seurat.R
# Load Xu et al. 2024 Primary Breast Tumor Atlas, subset to HR+ treatment-naive,
# apply QC, SCTransform, Harmony, PCA/UMAP
#
# Dataset: Xu et al. 2024 — 36 HR+ treatment-naive patients from 9 studies (~115K cells)
# Previously: Wu et al. (GSE176078) — 10 ER+ patients (~31K cells)
#
# Inputs:
#   - /ix1/alee/LO_LAB/General/Public_Data/BC-Datasets/Xu_etal_Primary_Breast_Tumor_Atlas_2024/matrix.mtx
#   - /ix1/alee/LO_LAB/General/Public_Data/BC-Datasets/Xu_etal_Primary_Breast_Tumor_Atlas_2024/genes.tsv
#   - /ix1/alee/LO_LAB/General/Public_Data/BC-Datasets/Xu_etal_Primary_Breast_Tumor_Atlas_2024/barcodes.tsv
#   - /ix1/alee/LO_LAB/General/Public_Data/BC-Datasets/Xu_etal_Primary_Breast_Tumor_Atlas_2024/metadata.csv
#
# Outputs:
#   - data/human_scrnaseq/SeuratObj_Xu2024_HRpos_AfterQCSCT.rds

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(glmGamPoi)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
PCA_DIMS <- 30

# QC thresholds (same as manuscript methods)
MIN_FEATURES <- 200
MAX_FEATURES <- 6000
MIN_COUNTS   <- 400
MAX_MITO_PCT <- 15

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
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
raw_data_dir <- "/ix1/alee/LO_LAB/General/Public_Data/BC-Datasets/Xu_etal_Primary_Breast_Tumor_Atlas_2024"
output_dir <- file.path(project_root, "data/human_scrnaseq")
log_dir <- file.path(script_dir, "outputs/preprocessing")

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Loading Xu et al. 2024 Primary Breast Tumor Atlas ===\n")
cat("Raw data directory:", raw_data_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# -----------------------------------------------------------------------------
# Step 1: Load sparse matrix and metadata
# -----------------------------------------------------------------------------
cat("Step 1: Loading sparse matrix...\n")

matrix_file <- file.path(raw_data_dir, "matrix.mtx")
genes_file <- file.path(raw_data_dir, "genes.tsv")
barcodes_file <- file.path(raw_data_dir, "barcodes.tsv")
metadata_file <- file.path(raw_data_dir, "metadata.csv")

counts <- readMM(matrix_file)

# Read gene names and deduplicate (58K genes may have duplicates)
# genes.tsv has two tab-separated columns (ID\tSymbol); take the symbol (2nd col)
genes_raw <- sapply(strsplit(readLines(genes_file), "\t"), function(x) x[length(x)])
genes <- make.unique(genes_raw)
cat("  Gene duplicates resolved:", sum(genes != genes_raw), "renamed\n")

# Use readLines (NOT read.table) — some barcodes contain spaces
# (e.g., "Aziziimmune_BC5_  260") which read.table splits into extra rows
barcodes <- readLines(barcodes_file)

rownames(counts) <- genes
colnames(counts) <- barcodes

cat("  Matrix dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")

# Load metadata — use read.csv (NOT fread) due to embedded commas in quoted fields
cat("  Loading metadata with read.csv (handles embedded commas)...\n")
metadata <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)

cat("  Metadata loaded:", nrow(metadata), "cells,", ncol(metadata), "columns\n\n")

# -----------------------------------------------------------------------------
# Step 2: Subset to HR+ treatment-naive patients
# -----------------------------------------------------------------------------
cat("Step 2: Subsetting to HR+ treatment-naive patients...\n")

# Filter by Clinical_Subtype == "HR+" and Treatment_Status == "Naive"
hr_mask <- metadata$Clinical_Subtype == "HR+" & metadata$Treatment_Status == "Naive"
cat("  HR+ Naive cells:", sum(hr_mask), "of", nrow(metadata), "total\n")

metadata_hr <- metadata[hr_mask, ]
counts_hr <- counts[, rownames(metadata_hr)]

cat("  HR+ patients:", length(unique(metadata_hr$Patient_ID)), "\n")
cat("  Datasets represented:", length(unique(metadata_hr$Dataset)), "\n")
cat("  Patient IDs:\n")
print(sort(table(metadata_hr$Patient_ID), decreasing = TRUE))

# -----------------------------------------------------------------------------
# Step 3: Create Seurat object
# -----------------------------------------------------------------------------
cat("\nStep 3: Creating Seurat object...\n")

# Set orig.ident to Patient_ID
# Do NOT use nCount_RNA/nFeature_RNA from metadata (fractional/pre-normalized)
# — let Seurat recompute from raw integer counts
meta_for_seurat <- metadata_hr[, !colnames(metadata_hr) %in% c("nCount_RNA", "nFeature_RNA")]

seurat_obj <- CreateSeuratObject(
  counts = counts_hr,
  meta.data = meta_for_seurat,
  project = "Xu2024_HRpos"
)

# Set orig.ident to Patient_ID
seurat_obj$orig.ident <- seurat_obj$Patient_ID

cat("  Seurat object created\n")
cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  Total genes:", nrow(seurat_obj), "\n\n")

# -----------------------------------------------------------------------------
# Step 4: QC filtering
# -----------------------------------------------------------------------------
cat("Step 4: QC filtering...\n")
cat("  Thresholds: nFeature [", MIN_FEATURES, ",", MAX_FEATURES, "],",
    "nCount >=", MIN_COUNTS, ", percent.mt <", MAX_MITO_PCT, "%\n")

# Calculate percent.mito from MT- genes (Xu doesn't provide pre-computed)
seurat_obj[["percent.mito"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

cat("  Before QC:", ncol(seurat_obj), "cells\n")

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= MIN_FEATURES &
           nFeature_RNA <= MAX_FEATURES &
           nCount_RNA >= MIN_COUNTS &
           percent.mito < MAX_MITO_PCT
)

cat("  After QC:", ncol(seurat_obj), "cells\n\n")

# Free memory from full matrix
rm(counts, counts_hr, metadata, metadata_hr, meta_for_seurat)
gc()

# -----------------------------------------------------------------------------
# Step 5: SCTransform normalization (with glmGamPoi for speed on 115K cells)
# -----------------------------------------------------------------------------
cat("Step 5: Running SCTransform (method = glmGamPoi)...\n")

# Increase future globals limit for large datasets (Xu atlas ~115K cells)
options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB

seurat_obj <- SCTransform(seurat_obj, method = "glmGamPoi", verbose = FALSE)

cat("  SCTransform complete\n\n")

# -----------------------------------------------------------------------------
# Step 6: PCA
# -----------------------------------------------------------------------------
cat("Step 6: Running PCA (", PCA_DIMS, " components)...\n")

seurat_obj <- RunPCA(seurat_obj, npcs = PCA_DIMS, verbose = FALSE)

# -----------------------------------------------------------------------------
# Step 7: Harmony integration (by Dataset — 9 studies need batch correction)
# -----------------------------------------------------------------------------
cat("Step 7: Running Harmony integration (group.by = Dataset)...\n")

seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "Dataset",
  assay.use = "SCT",
  plot_convergence = FALSE
)

cat("  Harmony complete\n\n")

# -----------------------------------------------------------------------------
# Step 8: UMAP (using Harmony-corrected embeddings)
# -----------------------------------------------------------------------------
cat("Step 8: Running UMAP on Harmony embeddings...\n")

seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "harmony",
  dims = 1:PCA_DIMS,
  verbose = FALSE
)

cat("  UMAP complete\n\n")

# -----------------------------------------------------------------------------
# Step 9: RNA assay NormalizeData + ScaleData
# Creates scale.data slot needed by downstream violin plots
# -----------------------------------------------------------------------------
cat("Step 9: NormalizeData + ScaleData on RNA assay...\n")

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

# Scale variable features + specific markers used with slot='scale.data' downstream
# (Full 59K x 115K scale.data would be ~50 GB — too expensive)
markers_to_scale <- unique(c(
  VariableFeatures(seurat_obj),
  # Macrophage M2 markers (05_gene_expression_violin.R, slot='scale.data')
  "CCL2", "TGFB1", "CD163", "MRC1",
  # Key markers for feature plots
  "EPCAM", "KRT19", "CD68", "CD3D", "MS4A1", "PECAM1",
  "CTLA4", "PDCD1"
))
markers_to_scale <- markers_to_scale[markers_to_scale %in% rownames(seurat_obj)]
cat("  Scaling", length(markers_to_scale), "features (variable + key markers)...\n")
seurat_obj <- ScaleData(seurat_obj, features = markers_to_scale, verbose = FALSE)

cat("  RNA assay: data and scale.data slots populated\n")
cat("  Normalized data:", nrow(GetAssayData(seurat_obj, assay = "RNA", slot = "data")), "genes x",
    ncol(GetAssayData(seurat_obj, assay = "RNA", slot = "data")), "cells\n\n")

# Reset default assay to SCT
DefaultAssay(seurat_obj) <- "SCT"

# -----------------------------------------------------------------------------
# Step 10: Map metadata columns for downstream compatibility
# -----------------------------------------------------------------------------
cat("Step 10: Mapping metadata columns...\n")

# Map Xu Cell_Type_Annotation → CellTypeMajor (used by downstream scripts)
seurat_obj$CellTypeMajor <- seurat_obj$Cell_Type_Annotation
seurat_obj$CellTypeMinor <- NA_character_
seurat_obj$Subtype <- seurat_obj$Molecular_Subtype

cat("  Cell type major:\n")
print(table(seurat_obj$CellTypeMajor))

# -----------------------------------------------------------------------------
# Step 11: Save output (BEFORE plots so a plot crash doesn't lose computation)
# -----------------------------------------------------------------------------
cat("\nStep 11: Saving processed object...\n")

output_file <- file.path(output_dir, "SeuratObj_Xu2024_HRpos_AfterQCSCT.rds")
saveRDS(seurat_obj, output_file)

cat("  Saved:", output_file, "\n")
cat("  Final dimensions:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")

# -----------------------------------------------------------------------------
# Step 12: Generate QC plots (non-critical, wrapped in tryCatch)
# -----------------------------------------------------------------------------
cat("\nStep 12: Generating QC plots...\n")

tryCatch({
  # UMAP by Dataset (post-Harmony) — check inter-study mixing
  p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "Dataset") +
    ggtitle("UMAP by Dataset (post-Harmony)")
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "CellTypeMajor", label = TRUE) +
    ggtitle("UMAP by Cell Type") +
    NoLegend()

  pdf(file.path(log_dir, "preprocessing_qc.pdf"), width = 14, height = 6)
  print(wrap_plots(list(p1, p2), ncol = 2))
  dev.off()

  # QC violin plots (group by Dataset instead of individual patients — too many)
  pdf(file.path(log_dir, "qc_violins.pdf"), width = 12, height = 4)
  qc_plots <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mito"), function(feat) {
    VlnPlot(seurat_obj, features = feat, group.by = "Dataset", pt.size = 0) +
      theme(legend.position = "none")
  })
  print(wrap_plots(qc_plots, ncol = 3))
  dev.off()

  cat("  Saved QC plots to:", log_dir, "\n")
}, error = function(e) {
  cat("  WARNING: QC plot generation failed:", e$message, "\n")
  cat("  (Non-critical — processed object already saved)\n")
  tryCatch(dev.off(), error = function(e2) NULL)
})

cat("\n=== Preprocessing complete ===\n")

# Print session info for reproducibility
cat("\nSession Info:\n")
print(sessionInfo())
