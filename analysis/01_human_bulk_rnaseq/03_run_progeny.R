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
# Step 1: Load TPM Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading TPM data...\n")
tpm_file <- file.path(data_dir, "raw/HumanERpAge_39404g168s_TPMlog2.txt")
gene_annot_file <- file.path(data_dir, "external/Homo_sapiens.gene_info.txt")
sample_annot_file <- file.path(output_dir, "sample_annotation.rds")

check_file_exists(tpm_file, "TPM data file")
check_file_exists(gene_annot_file, "Gene annotation file")
check_file_exists(sample_annot_file, "Sample annotation file")

tpm_data <- fread(tpm_file, stringsAsFactors = FALSE, header = TRUE)
colnames(tpm_data)[1] <- "GeneSymb"
colnames(tpm_data) <- gsub("_LEE.*", "", colnames(tpm_data))

# Load gene annotation for filtering
gene_annot <- fread(gene_annot_file, header = TRUE, stringsAsFactors = FALSE)
prot_coding <- gene_annot %>%
  filter(type_of_gene == "protein-coding", !grepl("^LOC\\d+", Symbol)) %>%
  pull(Symbol)

tpm_filtered <- tpm_data %>%
  filter(GeneSymb %in% prot_coding)

cat("  TPM matrix:", nrow(tpm_filtered), "genes x", ncol(tpm_filtered) - 1, "samples\n")

# -----------------------------------------------------------------------------
# Step 2: Prepare Expression Matrix
# -----------------------------------------------------------------------------
cat("Step 2: Preparing expression matrix...\n")
sample_annot <- readRDS(sample_annot_file)

tpm_t <- tpm_filtered %>%
  column_to_rownames("GeneSymb") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleName")

tpm_annot <- inner_join(
  sample_annot[, c("SampleName", "SampleNameGroup")],
  tpm_t,
  by = "SampleName"
) %>%
  select(-SampleName) %>%
  column_to_rownames("SampleNameGroup") %>%
  t() %>%
  as.data.frame()

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
