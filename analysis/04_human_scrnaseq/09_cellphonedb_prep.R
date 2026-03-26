#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/09_cellphonedb_prep.R
# Prepare input files for CellPhoneDB analysis
#
# Gene filtering matches Sanghoon's original code:
#   Young/Elderly: rowSums > 20 (CellPhoneDB-specific)
#   Whole cohort:  rowSums > 10, remove non-protein-coding (contains "."),
#                  colSums > 1000 → ~18,063 genes x ~28,732 cells
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/cellphonedb/counts_*.txt
#   - analysis/04_human_scrnaseq/outputs/cellphonedb/metadata_*.txt

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
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
output_dir <- file.path(script_dir, "outputs/cellphonedb")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== CellPhoneDB Input Preparation ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_file <- file.path(dirname(output_dir), "seurat_annotated.rds")
seurat_obj <- readRDS(seurat_file)

cat("  Cells:", ncol(seurat_obj), "\n")

# Use CellTypeAnnot as identity (matching 03_cell_type_annotation.R)
Idents(seurat_obj) <- seurat_obj$CellTypeAnnot

seurat_young <- subset(seurat_obj, AgeGroup == "Young")
seurat_elderly <- subset(seurat_obj, AgeGroup == "Elderly")
seurat_whole <- seurat_obj

# -----------------------------------------------------------------------------
# Step 2: Functions to prepare CellPhoneDB inputs
# -----------------------------------------------------------------------------

# Young/Elderly: simple rowSums > 20 filter (matches original)
prepare_cellphonedb_agegroup <- function(seurat_subset, name) {
  cat(sprintf("\nProcessing %s (%d cells)...\n", name, ncol(seurat_subset)))

  counts <- GetAssayData(seurat_subset, slot = "counts")

  # Gene filter: rowSums > 20 (matches original Young/Elderly)
  counts_filtered <- counts[rowSums(counts) > 20, ]
  cat(sprintf("  After gene filter (rowSums > 20): %d genes x %d cells\n",
              nrow(counts_filtered), ncol(counts_filtered)))

  # Replace hyphens in cell IDs (CellPhoneDB requirement)
  colnames(counts_filtered) <- gsub("-", "_", colnames(counts_filtered))

  # Save counts
  count_file <- file.path(output_dir, sprintf("counts_%s.txt", name))
  fwrite(
    as.data.frame(as.matrix(counts_filtered)) %>% tibble::rownames_to_column("Gene"),
    count_file, sep = "\t", quote = FALSE
  )
  cat(sprintf("  Saved: %s\n", basename(count_file)))

  # Metadata
  cells_keep <- colnames(counts_filtered)
  original_cells <- gsub("_", "-", cells_keep)
  matched_idx <- match(original_cells, colnames(seurat_subset))

  metadata <- data.frame(
    Cell = cells_keep,
    cell_type = gsub("-|_|\\+", "", seurat_subset$CellTypeAnnot[matched_idx])
  )

  meta_file <- file.path(output_dir, sprintf("metadata_%s.txt", name))
  fwrite(metadata, meta_file, sep = "\t", quote = FALSE)
  cat(sprintf("  Saved: %s\n", basename(meta_file)))

  return(list(n_genes = nrow(counts_filtered), n_cells = ncol(counts_filtered)))
}

# Whole cohort: matches original's 3-step filtering
# Original code:
#   1. rowSums(counts) > 10
#   2. Remove non-protein-coding genes (containing "." in gene name)
#   3. colSums > 1000
#   → 18,063 genes x 28,732 cells
prepare_cellphonedb_whole <- function(seurat_subset, name) {
  cat(sprintf("\nProcessing %s (%d cells)...\n", name, ncol(seurat_subset)))

  counts <- GetAssayData(seurat_subset, slot = "counts")
  counts_df <- as.data.frame(as.matrix(counts))

  cat(sprintf("  Starting: %d genes x %d cells\n", nrow(counts_df), ncol(counts_df)))

  # Step 1: Gene filter — rowSums > 10 (NOT 20, matches original whole-cohort code)
  counts_df <- counts_df[rowSums(counts_df) > 10, ]
  cat(sprintf("  After gene filter (rowSums > 10): %d genes\n", nrow(counts_df)))

  # Step 2: Remove non-protein-coding genes (containing "." in gene name)
  # Original: dplyr::filter(!grepl("\\.", rownames(...)))
  non_coding <- grepl("\\.", rownames(counts_df))
  counts_df <- counts_df[!non_coding, ]
  cat(sprintf("  After removing non-protein-coding: %d genes\n", nrow(counts_df)))

  # Step 3: Cell filter — colSums > 1000
  counts_df <- counts_df[, colSums(counts_df) > 1000]
  cat(sprintf("  After cell filter (colSums > 1000): %d genes x %d cells\n",
              nrow(counts_df), ncol(counts_df)))
  cat("  (Original expects: ~18,063 genes x ~28,732 cells)\n")

  # Replace hyphens
  colnames(counts_df) <- gsub("-", "_", colnames(counts_df))

  # Save counts
  count_file <- file.path(output_dir, sprintf("counts_%s.txt", name))
  fwrite(
    counts_df %>% tibble::rownames_to_column("Gene"),
    count_file, sep = "\t", quote = FALSE
  )
  cat(sprintf("  Saved: %s\n", basename(count_file)))

  # Metadata — subset to surviving cells
  cells_keep <- colnames(counts_df)
  original_cells <- gsub("_", "-", cells_keep)
  matched_idx <- match(original_cells, colnames(seurat_subset))

  metadata <- data.frame(
    Cell = cells_keep,
    cell_type = gsub("-|_|\\+", "", seurat_subset$CellTypeAnnot[matched_idx])
  )

  meta_file <- file.path(output_dir, sprintf("metadata_%s.txt", name))
  fwrite(metadata, meta_file, sep = "\t", quote = FALSE)
  cat(sprintf("  Saved: %s\n", basename(meta_file)))

  return(list(n_genes = nrow(counts_df), n_cells = ncol(counts_df)))
}

# -----------------------------------------------------------------------------
# Step 3: Prepare files for each group
# -----------------------------------------------------------------------------
cat("\nStep 3: Preparing CellPhoneDB inputs...\n")

stats_young <- prepare_cellphonedb_agegroup(seurat_young, "Young")
stats_elderly <- prepare_cellphonedb_agegroup(seurat_elderly, "Elderly")
stats_whole <- prepare_cellphonedb_whole(seurat_whole, "Whole")

# -----------------------------------------------------------------------------
# Step 4: Summary
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat(sprintf("Young:   %d genes x %d cells (filter: rowSums > 20)\n",
            stats_young$n_genes, stats_young$n_cells))
cat(sprintf("Elderly: %d genes x %d cells (filter: rowSums > 20)\n",
            stats_elderly$n_genes, stats_elderly$n_cells))
cat(sprintf("Whole:   %d genes x %d cells (filter: rowSums > 10, no '.', colSums > 1000)\n",
            stats_whole$n_genes, stats_whole$n_cells))

cat("\n=== CellPhoneDB prep complete ===\n")
cat("Run CellPhoneDB with:\n")
cat("  cellphonedb method statistical_analysis metadata_<group>.txt counts_<group>.txt\n")
