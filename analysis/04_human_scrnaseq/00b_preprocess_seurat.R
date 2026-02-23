#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/00b_preprocess_seurat.R
# Preprocess Wu et al. (GSE176078) scRNA-seq data following published methods
#
# Methods citation:
#   QC: nFeature 200-6000, nCount >400, percent.mt <15%
#   Normalization: SCTransform per sample (glmGamPoi)
#   Dimensionality reduction: PCA (30 dims) -> UMAP
#   Cell types: Original Wu et al. annotations
#
# Inputs:
#   - data/human_scrnaseq/raw/CID*/matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
#
# Outputs:
#   - data/human_scrnaseq/SeuratObj_GSE176078_ERpos_AfterQCSCT.rds

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(future)
})

# Enable parallel processing for SCTransform
plan("multicore", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

# -----------------------------------------------------------------------------
# Configuration (from published methods)
# -----------------------------------------------------------------------------
MIN_FEATURES <- 200
MAX_FEATURES <- 6000
MIN_COUNTS <- 400
MAX_MT_PERCENT <- 15
PCA_DIMS <- 30

# Sample IDs (10 ER+ samples)
SAMPLE_IDS <- c(
  "CID3941", "CID4530N", "CID4535",           # Young
  "CID4463", "CID4040", "CID4471", "CID4461", # MidAge
  "CID3948", "CID4067", "CID4290A"            # Elderly
)

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
raw_data_dir <- file.path(project_root, "data/human_scrnaseq/raw")
output_dir <- file.path(project_root, "data/human_scrnaseq")
log_dir <- file.path(script_dir, "outputs/preprocessing")

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Wu et al. scRNA-seq Preprocessing ===\n")
cat("Raw data directory:", raw_data_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# -----------------------------------------------------------------------------
# Step 1: Load and QC each sample
# -----------------------------------------------------------------------------
cat("Step 1: Loading and QC filtering samples...\n")

qc_stats <- data.frame(
  sample = character(),
  cells_raw = integer(),
  cells_after_qc = integer(),
  stringsAsFactors = FALSE
)

seurat_list <- lapply(SAMPLE_IDS, function(sample_id) {
  cat("  Processing:", sample_id, "... ")

  sample_dir <- file.path(raw_data_dir, sample_id)

  if (!dir.exists(sample_dir)) {
    stop(paste("Sample directory not found:", sample_dir))
  }

  # Load 10X data
  counts <- Read10X(data.dir = sample_dir)
  obj <- CreateSeuratObject(counts = counts, project = sample_id)
  obj$orig.ident <- sample_id

  cells_raw <- ncol(obj)

  # QC metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  # QC filtering (per methods)
  obj <- subset(obj,
    subset = nFeature_RNA > MIN_FEATURES &
             nFeature_RNA < MAX_FEATURES &
             nCount_RNA > MIN_COUNTS &
             percent.mt < MAX_MT_PERCENT
  )

  cells_after_qc <- ncol(obj)

  cat(cells_raw, "->", cells_after_qc, "cells\n")

  # Track QC stats
  qc_stats <<- rbind(qc_stats, data.frame(
    sample = sample_id,
    cells_raw = cells_raw,
    cells_after_qc = cells_after_qc
  ))

  # SCTransform per sample (for batch correction)
  obj <- SCTransform(obj,
    method = "glmGamPoi",
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )

  return(obj)
})

names(seurat_list) <- SAMPLE_IDS

# Save QC stats
fwrite(qc_stats, file.path(log_dir, "qc_stats.tsv"), sep = "\t")
cat("\nQC Summary:\n")
print(qc_stats)
cat("Total cells after QC:", sum(qc_stats$cells_after_qc), "\n\n")

# -----------------------------------------------------------------------------
# Step 2: Merge samples
# -----------------------------------------------------------------------------
cat("Step 2: Merging samples...\n")

merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = SAMPLE_IDS,
  project = "Wu_GSE176078_ERpos"
)

cat("  Merged object:", ncol(merged), "cells x", nrow(merged), "genes\n\n")

# -----------------------------------------------------------------------------
# Step 3: PCA and UMAP
# -----------------------------------------------------------------------------
cat("Step 3: Running PCA and UMAP...\n")

# Need to re-run variable features on merged object
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = PCA_DIMS, verbose = FALSE)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:PCA_DIMS, verbose = FALSE)

cat("  PCA dims:", PCA_DIMS, "\n")
cat("  UMAP computed\n\n")

# -----------------------------------------------------------------------------
# Step 4: Generate QC plots
# -----------------------------------------------------------------------------
cat("Step 4: Generating QC plots...\n")

# UMAP by sample
p1 <- DimPlot(merged, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP by Sample")

# UMAP density
p2 <- DimPlot(merged, reduction = "umap") +
  ggtitle("UMAP All Cells")

pdf(file.path(log_dir, "preprocessing_qc.pdf"), width = 12, height = 5)
print(p1 + p2)
dev.off()

# Violin plots of QC metrics
pdf(file.path(log_dir, "qc_violins.pdf"), width = 12, height = 4)
print(VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "orig.ident", ncol = 3, pt.size = 0))
dev.off()

cat("  Saved QC plots to:", log_dir, "\n\n")

# -----------------------------------------------------------------------------
# Step 5: Save output
# -----------------------------------------------------------------------------
cat("Step 5: Saving processed object...\n")

output_file <- file.path(output_dir, "SeuratObj_GSE176078_ERpos_AfterQCSCT.rds")
saveRDS(merged, output_file)

cat("  Saved:", output_file, "\n")
cat("  Final dimensions:", ncol(merged), "cells x", nrow(merged), "genes\n")

cat("\n=== Preprocessing complete ===\n")

# Print session info for reproducibility
cat("\nSession Info:\n")
print(sessionInfo())
