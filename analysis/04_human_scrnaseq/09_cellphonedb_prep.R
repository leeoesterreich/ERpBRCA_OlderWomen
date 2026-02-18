#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/09_cellphonedb_prep.R
# Prepare input files for CellPhoneDB analysis
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

# FIX: Define these subsets properly instead of using undefined variables
seurat_young <- subset(seurat_obj, AgeGroup == "Young")
seurat_elderly <- subset(seurat_obj, AgeGroup == "Elderly")
seurat_whole <- seurat_obj

# -----------------------------------------------------------------------------
# Step 2: Function to prepare CellPhoneDB inputs
# -----------------------------------------------------------------------------
prepare_cellphonedb <- function(seurat_subset, name) {
  cat(sprintf("\nProcessing %s (%d cells)...\n", name, ncol(seurat_subset)))

  # Get count data
  counts <- GetAssayData(seurat_subset, slot = "counts")

  # Filter low-count genes (>20 total counts across cells)
  counts_filtered <- counts[rowSums(counts) > 20, ]

  # Filter low-count cells (>1000 total counts)
  counts_filtered <- counts_filtered[, colSums(counts_filtered) > 1000]

  # FIX: Replace hyphens in cell IDs (CellPhoneDB requirement)
  colnames(counts_filtered) <- gsub("-", "_", colnames(counts_filtered))

  # Save counts
  count_file <- file.path(output_dir, sprintf("counts_%s.txt", name))
  fwrite(
    as.data.frame(counts_filtered) %>% tibble::rownames_to_column("Gene"),
    count_file,
    sep = "\t",
    quote = FALSE
  )
  cat(sprintf("  Saved: %s (%d genes x %d cells)\n",
              basename(count_file), nrow(counts_filtered), ncol(counts_filtered)))

  # Prepare metadata
  cells_keep <- colnames(counts_filtered)
  # FIX: Match cell names after hyphen replacement
  original_cells <- gsub("_", "-", cells_keep)

  metadata <- data.frame(
    Cell = cells_keep,
    cell_type = gsub("-|_|\\+", "", seurat_subset$CellTypeAnnot[match(original_cells, colnames(seurat_subset))])
  )

  # FIX: Verify cell IDs match between counts and metadata
  if (!all(metadata$Cell == colnames(counts_filtered))) {
    warning("Cell ID mismatch between counts and metadata!")
  }

  # Save metadata
  meta_file <- file.path(output_dir, sprintf("metadata_%s.txt", name))
  fwrite(metadata, meta_file, sep = "\t", quote = FALSE)
  cat(sprintf("  Saved: %s\n", basename(meta_file)))

  return(list(
    n_genes = nrow(counts_filtered),
    n_cells = ncol(counts_filtered)
  ))
}

# -----------------------------------------------------------------------------
# Step 3: Prepare files for each group
# -----------------------------------------------------------------------------
cat("\nStep 3: Preparing CellPhoneDB inputs...\n")

stats_young <- prepare_cellphonedb(seurat_young, "Young")
stats_elderly <- prepare_cellphonedb(seurat_elderly, "Elderly")
stats_whole <- prepare_cellphonedb(seurat_whole, "Whole")

# -----------------------------------------------------------------------------
# Step 4: Summary
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat(sprintf("Young: %d genes x %d cells\n", stats_young$n_genes, stats_young$n_cells))
cat(sprintf("Elderly: %d genes x %d cells\n", stats_elderly$n_genes, stats_elderly$n_cells))
cat(sprintf("Whole: %d genes x %d cells\n", stats_whole$n_genes, stats_whole$n_cells))

cat("\n=== CellPhoneDB prep complete ===\n")
cat("Run CellPhoneDB with:\n")
cat("  cellphonedb method statistical_analysis metadata_<group>.txt counts_<group>.txt\n")
