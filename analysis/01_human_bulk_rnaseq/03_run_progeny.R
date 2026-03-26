#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/03_run_progeny.R
# Run PROGENy pathway activity analysis
#
# Inputs:
#   - data/human_bulk_rnaseq/raw/HumanERpAge_39404g168s_TPMlog2.txt
#   - analysis/01_human_bulk_rnaseq/outputs/sample_annotation.rds
#
# Outputs:
#   - analysis/01_human_bulk_rnaseq/outputs/progeny_pathway_activity.rds

# BIOSTATISTICAL FIX: Set random seed for reproducibility
set.seed(12345)

suppressPackageStartupMessages({
  library(progeny)
  library(dplyr)
  library(data.table)
  library(tibble)
})

# Helper function to check file existence
check_file_exists <- function(filepath, description = "file") {
  if (!file.exists(filepath)) {
    stop(sprintf("ERROR: %s not found: %s", description, filepath))
  }
  cat(sprintf("  Found: %s\n", basename(filepath)))
}

# Define paths - use commandArgs to get script directory when run via Rscript
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
data_dir <- file.path(project_root, "data/human_bulk_rnaseq")

cat("=== PROGENy Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Preprocessed Expression Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading preprocessed expression data...\n")
vst_file <- file.path(output_dir, "vst_normalized_matrix.rds")
sample_annot_file <- file.path(output_dir, "sample_annotation.rds")

check_file_exists(vst_file, "VST-normalized matrix (from 01_preprocess.R)")
check_file_exists(sample_annot_file, "Sample annotation file")

# Use the same preprocessed matrix as GSVA (consistency)
vst_matrix <- readRDS(vst_file)
sample_annot <- readRDS(sample_annot_file)

cat("  VST matrix:", nrow(vst_matrix), "genes x", ncol(vst_matrix), "samples\n")

# -----------------------------------------------------------------------------
# Step 2: Prepare Expression Matrix
# -----------------------------------------------------------------------------
cat("Step 2: Preparing expression matrix...\n")

# Convert to base data.frame if data.table (preserves row names for PROGENy)
if (inherits(vst_matrix, "data.table")) {
  gene_names <- vst_matrix[[1]]  # first column is typically gene names
  # Check if first column is character (gene IDs)
  if (is.character(gene_names)) {
    vst_matrix <- as.data.frame(vst_matrix[, -1])
    rownames(vst_matrix) <- gene_names
  } else {
    vst_matrix <- as.data.frame(vst_matrix)
  }
}

# Align sample names with annotation
common_samples <- intersect(colnames(vst_matrix), sample_annot$SampleName)
if (length(common_samples) == 0) {
  # Try SampleNameGroup matching
  common_samples <- intersect(colnames(vst_matrix), sample_annot$SampleNameGroup)
  tpm_annot <- vst_matrix[, common_samples, drop = FALSE]
} else {
  # Rename to SampleNameGroup for downstream compatibility
  name_map <- setNames(sample_annot$SampleNameGroup, sample_annot$SampleName)
  tpm_annot <- vst_matrix[, common_samples, drop = FALSE]
  colnames(tpm_annot) <- name_map[common_samples]
}

cat("  Expression matrix:", nrow(tpm_annot), "genes x", ncol(tpm_annot), "samples\n")

# -----------------------------------------------------------------------------
# Step 3: Run PROGENy
# -----------------------------------------------------------------------------
cat("Step 3: Running PROGENy...\n")

# BIOSTATISTICAL FIX: Random seed already set at script start
pathway_activity <- progeny(
  as.matrix(tpm_annot),
  scale = FALSE,
  organism = "Human",
  top = 100,
  perm = 1000
)

cat("  Pathway activity:", nrow(pathway_activity), "samples x", ncol(pathway_activity), "pathways\n")

# -----------------------------------------------------------------------------
# Step 4: Save Outputs
# -----------------------------------------------------------------------------
cat("Step 4: Saving outputs...\n")
saveRDS(pathway_activity, file.path(output_dir, "progeny_pathway_activity.rds"))

cat("\n=== PROGENy complete ===\n")
