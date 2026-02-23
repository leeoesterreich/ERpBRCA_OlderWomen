#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/00b_preprocess_seurat.R
# Load Wu et al. (GSE176078) pre-processed scRNA-seq data
#
# The GEO data is provided as a combined sparse matrix with metadata
# including QC metrics and cell type annotations from Wu et al.
#
# Inputs:
#   - data/human_scrnaseq/raw/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx
#   - data/human_scrnaseq/raw/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv
#   - data/human_scrnaseq/raw/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv
#   - data/human_scrnaseq/raw/Wu_etal_2021_BRCA_scRNASeq/metadata.csv
#
# Outputs:
#   - data/human_scrnaseq/SeuratObj_GSE176078_ERpos_AfterQCSCT.rds

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
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
raw_data_dir <- file.path(project_root, "data/human_scrnaseq/raw/Wu_etal_2021_BRCA_scRNASeq")
output_dir <- file.path(project_root, "data/human_scrnaseq")
log_dir <- file.path(script_dir, "outputs/preprocessing")

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Loading Wu et al. Pre-processed scRNA-seq Data ===\n")
cat("Raw data directory:", raw_data_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# -----------------------------------------------------------------------------
# Step 1: Load sparse matrix and metadata
# -----------------------------------------------------------------------------
cat("Step 1: Loading sparse matrix...\n")

# Load matrix components
matrix_file <- file.path(raw_data_dir, "count_matrix_sparse.mtx")
genes_file <- file.path(raw_data_dir, "count_matrix_genes.tsv")
barcodes_file <- file.path(raw_data_dir, "count_matrix_barcodes.tsv")
metadata_file <- file.path(raw_data_dir, "metadata.csv")

counts <- readMM(matrix_file)
genes <- fread(genes_file, header = FALSE)$V1
barcodes <- fread(barcodes_file, header = FALSE)$V1

# Set matrix dimensions
rownames(counts) <- genes
colnames(counts) <- barcodes

cat("  Matrix dimensions:", nrow(counts), "genes x", ncol(counts), "cells\n")

# Load metadata
metadata <- fread(metadata_file, header = TRUE)
rownames(metadata) <- metadata$V1
metadata$V1 <- NULL

cat("  Metadata loaded:", nrow(metadata), "cells\n\n")

# -----------------------------------------------------------------------------
# Step 2: Create Seurat object
# -----------------------------------------------------------------------------
cat("Step 2: Creating Seurat object...\n")

seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = as.data.frame(metadata),
  project = "Wu_GSE176078"
)

cat("  Seurat object created\n")
cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  Total genes:", nrow(seurat_obj), "\n\n")

# -----------------------------------------------------------------------------
# Step 3: Subset to 10 ER+ samples
# -----------------------------------------------------------------------------
cat("Step 3: Subsetting to 10 ER+ samples...\n")

seurat_subset <- subset(seurat_obj, subset = orig.ident %in% SAMPLE_IDS)

cat("  Cells after subset:", ncol(seurat_subset), "\n")
cat("  Samples:\n")
print(table(seurat_subset$orig.ident))

# -----------------------------------------------------------------------------
# Step 4: Normalize and run dimensionality reduction
# -----------------------------------------------------------------------------
cat("\nStep 4: Running SCTransform and dimensionality reduction...\n")

# SCTransform (data is already QC'd by Wu et al.)
seurat_subset <- SCTransform(seurat_subset,
                              method = "glmGamPoi",
                              vars.to.regress = "percent.mito",
                              verbose = FALSE)

# PCA
seurat_subset <- RunPCA(seurat_subset, npcs = PCA_DIMS, verbose = FALSE)

# UMAP
seurat_subset <- RunUMAP(seurat_subset,
                          reduction = "pca",
                          dims = 1:PCA_DIMS,
                          verbose = FALSE)

cat("  SCTransform, PCA, UMAP complete\n\n")

# -----------------------------------------------------------------------------
# Step 5: Rename cell type columns for consistency
# -----------------------------------------------------------------------------
cat("Step 5: Standardizing cell type annotations...\n")

# Wu et al. annotations
seurat_subset$CellTypeMajor <- seurat_subset$celltype_major
seurat_subset$CellTypeMinor <- seurat_subset$celltype_minor
seurat_subset$CellTypeSubset <- seurat_subset$celltype_subset
seurat_subset$Subtype <- seurat_subset$subtype

cat("  Cell type major:\n")
print(table(seurat_subset$CellTypeMajor))

# -----------------------------------------------------------------------------
# Step 6: Generate QC plots
# -----------------------------------------------------------------------------
cat("\nStep 6: Generating QC plots...\n")

# UMAP by sample
p1 <- DimPlot(seurat_subset, reduction = "umap", group.by = "orig.ident") +
  ggtitle("UMAP by Sample")

# UMAP by cell type
p2 <- DimPlot(seurat_subset, reduction = "umap", group.by = "CellTypeMajor", label = TRUE) +
  ggtitle("UMAP by Cell Type") +
  NoLegend()

pdf(file.path(log_dir, "preprocessing_qc.pdf"), width = 14, height = 6)
print(p1 + p2)
dev.off()

# Violin plots of QC metrics
pdf(file.path(log_dir, "qc_violins.pdf"), width = 12, height = 4)
print(VlnPlot(seurat_subset,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
              group.by = "orig.ident", ncol = 3, pt.size = 0))
dev.off()

cat("  Saved QC plots to:", log_dir, "\n\n")

# -----------------------------------------------------------------------------
# Step 7: Save output
# -----------------------------------------------------------------------------
cat("Step 7: Saving processed object...\n")

output_file <- file.path(output_dir, "SeuratObj_GSE176078_ERpos_AfterQCSCT.rds")
saveRDS(seurat_subset, output_file)

cat("  Saved:", output_file, "\n")
cat("  Final dimensions:", ncol(seurat_subset), "cells x", nrow(seurat_subset), "genes\n")

cat("\n=== Preprocessing complete ===\n")

# Print session info for reproducibility
cat("\nSession Info:\n")
print(sessionInfo())
