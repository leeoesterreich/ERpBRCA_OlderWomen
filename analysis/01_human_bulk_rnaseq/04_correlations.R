#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/04_correlations.R
# Calculate correlations between gene expression and pathway activity
# BIOSTATISTICAL FIX: Implements Benjamini-Hochberg FDR correction
#
# Inputs:
#   - analysis/01_human_bulk_rnaseq/outputs/progeny_pathway_activity.rds
#   - analysis/01_human_bulk_rnaseq/outputs/gsva_estrogen_pathways.rds
#   - data/human_bulk_rnaseq/raw/HumanERpAge_39404g168s_TPMlog2.txt
#
# Outputs:
#   - results/original/human_bulk_rnaseq/correlations_no_fdr.csv
#   - results/corrected/human_bulk_rnaseq/correlations_with_fdr.csv

set.seed(12345)

suppressPackageStartupMessages({
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

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
data_dir <- file.path(project_root, "data/human_bulk_rnaseq")
results_original <- file.path(project_root, "results/original/human_bulk_rnaseq")
results_corrected <- file.path(project_root, "results/corrected/human_bulk_rnaseq")

dir.create(results_original, showWarnings = FALSE, recursive = TRUE)
dir.create(results_corrected, showWarnings = FALSE, recursive = TRUE)

cat("=== Correlation Analysis with FDR Correction ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

# Define input files
progeny_file <- file.path(output_dir, "progeny_pathway_activity.rds")
gsva_file <- file.path(output_dir, "gsva_estrogen_pathways.rds")
sample_file <- file.path(output_dir, "sample_annotation.rds")
tpm_file <- file.path(data_dir, "raw/HumanERpAge_39404g168s_TPMlog2.txt")

# Check all input files exist
check_file_exists(progeny_file, "PROGENy pathway activity")
check_file_exists(gsva_file, "GSVA results")
check_file_exists(sample_file, "Sample annotation")
check_file_exists(tpm_file, "TPM data")

pathway_activity <- readRDS(progeny_file)
gsva_result <- readRDS(gsva_file)
sample_annot <- readRDS(sample_file)

# Load TPM data for gene expression
tpm_data <- fread(tpm_file, stringsAsFactors = FALSE, header = TRUE)
colnames(tpm_data)[1] <- "GeneSymb"
colnames(tpm_data) <- gsub("_LEE.*", "", colnames(tpm_data))

# -----------------------------------------------------------------------------
# Step 2: Define Genes and Pathways of Interest
# -----------------------------------------------------------------------------
cat("Step 2: Setting up analysis...\n")

# Genes of interest (from original analysis)
genes_of_interest <- c("PAK4", "HSD17B7", "GREB1", "PGR", "ESR1",
                       "TFF1", "CYP19A1", "HSD17B2", "SAA1", "Age")

# Filter to tumor samples only, excluding middle age
tumor_samples <- sample_annot %>%
  filter(Group == "Tumor", AgeRange %in% c("Young", "Elderly")) %>%
  pull(SampleNameGroup)

cat("  Genes:", length(genes_of_interest) - 1, "+ Age\n")
cat("  Tumor samples (Young + Elderly):", length(tumor_samples), "\n")

# -----------------------------------------------------------------------------
# Step 3: Prepare Combined Data
# -----------------------------------------------------------------------------
cat("Step 3: Preparing combined data...\n")

# Get Estrogen pathway activity from PROGENy
progeny_estrogen <- pathway_activity[tumor_samples, "Estrogen", drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("SampleNameGroup")

# Get GSVA scores
gsva_t <- gsva_result %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleNameGroup") %>%
  filter(SampleNameGroup %in% tumor_samples)

# Get gene expression
tpm_t <- tpm_data %>%
  column_to_rownames("GeneSymb") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleName")

gene_expr <- inner_join(
  sample_annot[, c("SampleName", "SampleNameGroup", "Age")],
  tpm_t,
  by = "SampleName"
) %>%
  filter(SampleNameGroup %in% tumor_samples) %>%
  select(SampleNameGroup, Age, any_of(genes_of_interest[genes_of_interest != "Age"]))

# Combine all data
combined_data <- gene_expr %>%
  inner_join(progeny_estrogen, by = "SampleNameGroup") %>%
  inner_join(gsva_t, by = "SampleNameGroup")

pathway_names <- c("Estrogen", colnames(gsva_t)[-1])
cat("  Combined data:", nrow(combined_data), "samples\n")
cat("  Pathways:", length(pathway_names), "\n")

# -----------------------------------------------------------------------------
# Step 4: Calculate Correlations
# -----------------------------------------------------------------------------
cat("Step 4: Calculating correlations...\n")

correlation_results <- data.frame()

for (gene in genes_of_interest) {
  for (pathway in pathway_names) {

    if (gene == "Age") {
      x_values <- combined_data$Age
    } else {
      if (!gene %in% colnames(combined_data)) next
      x_values <- combined_data[[gene]]
    }

    y_values <- combined_data[[pathway]]

    # Skip if insufficient data
    if (sum(!is.na(x_values) & !is.na(y_values)) < 5) next

    # Spearman correlation
    cor_test <- cor.test(x_values, y_values, method = "spearman", exact = FALSE)

    correlation_results <- rbind(correlation_results, data.frame(
      GeneSymb = gene,
      PathwayName = pathway,
      Spearman_Rho = cor_test$estimate,
      Spearman_pval = cor_test$p.value,
      n_samples = sum(!is.na(x_values) & !is.na(y_values))
    ))
  }
}

cat("  Total correlation tests:", nrow(correlation_results), "\n")

# -----------------------------------------------------------------------------
# Step 5: Apply FDR Correction
# -----------------------------------------------------------------------------
cat("Step 5: Applying FDR correction...\n")

# BIOSTATISTICAL FIX: Benjamini-Hochberg FDR correction
correlation_results$FDR_qval <- p.adjust(correlation_results$Spearman_pval, method = "BH")

# Add significance flags
correlation_results$Sig_nominal <- correlation_results$Spearman_pval < 0.05
correlation_results$Sig_FDR <- correlation_results$FDR_qval < 0.05

# Summary statistics
n_sig_nominal <- sum(correlation_results$Sig_nominal)
n_sig_fdr <- sum(correlation_results$Sig_FDR)
cat("  Significant at p < 0.05:", n_sig_nominal, "\n")
cat("  Significant at FDR < 0.05:", n_sig_fdr, "\n")

# -----------------------------------------------------------------------------
# Step 6: Save Results
# -----------------------------------------------------------------------------
cat("Step 6: Saving results...\n")

# Original (no FDR) - for comparison
original_results <- correlation_results %>%
  select(GeneSymb, PathwayName, Spearman_Rho, Spearman_pval, n_samples, Sig_nominal)

# Corrected (with FDR)
corrected_results <- correlation_results %>%
  select(GeneSymb, PathwayName, Spearman_Rho, Spearman_pval, FDR_qval,
         n_samples, Sig_nominal, Sig_FDR)

write.csv(original_results, file.path(results_original, "correlations_no_fdr.csv"),
          row.names = FALSE)
write.csv(corrected_results, file.path(results_corrected, "correlations_with_fdr.csv"),
          row.names = FALSE)

# Also save to analysis outputs
saveRDS(correlation_results, file.path(output_dir, "correlation_results.rds"))

cat("\n=== Correlation analysis complete ===\n")
cat("Original results:", file.path(results_original, "correlations_no_fdr.csv"), "\n")
cat("Corrected results:", file.path(results_corrected, "correlations_with_fdr.csv"), "\n")
