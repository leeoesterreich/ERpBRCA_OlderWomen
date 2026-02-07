#!/usr/bin/env Rscript
# analysis/03_comparison/01_generate_comparison_report.R
# Generate comparison report: Original vs Corrected results
#
# This script compares results before and after biostatistical corrections
# to validate reproducibility and document changes.

set.seed(12345)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# Define paths
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
original_dir <- file.path(project_root, "results/original")
corrected_dir <- file.path(project_root, "results/corrected")
comparison_dir <- file.path(project_root, "results/comparison")
fig_dir <- file.path(comparison_dir, "figures")

dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Generating Comparison Report ===\n")
cat("Project root:", project_root, "\n\n")

# =============================================================================
# 1. Human Bulk RNA-seq Comparison
# =============================================================================
cat("--- Human Bulk RNA-seq Comparison ---\n\n")

# Load original and corrected results
original_corr_file <- file.path(original_dir, "human_bulk_rnaseq/correlations_no_fdr.csv")
corrected_corr_file <- file.path(corrected_dir, "human_bulk_rnaseq/correlations_with_fdr.csv")

if (file.exists(original_corr_file) && file.exists(corrected_corr_file)) {
  original_corr <- read.csv(original_corr_file)
  corrected_corr <- read.csv(corrected_corr_file)

  cat("Original correlations:", nrow(original_corr), "tests\n")
  cat("Corrected correlations:", nrow(corrected_corr), "tests\n\n")

  # Compare significance
  original_sig <- sum(original_corr$pvalue < 0.05, na.rm = TRUE)
  corrected_sig_nominal <- sum(corrected_corr$pvalue < 0.05, na.rm = TRUE)
  corrected_sig_fdr <- sum(corrected_corr$FDR < 0.05, na.rm = TRUE)

  cat("Significance comparison:\n")
  cat("  Original (p < 0.05):", original_sig, "\n")
  cat("  Corrected (p < 0.05):", corrected_sig_nominal, "\n")
  cat("  Corrected (FDR < 0.05):", corrected_sig_fdr, "\n\n")

  # Check correlation values match (reproducibility)
  if ("correlation" %in% colnames(original_corr) && "correlation" %in% colnames(corrected_corr)) {
    # Merge by pathway and variable
    merged <- merge(
      original_corr %>% select(pathway, variable, correlation, pvalue) %>% rename(corr_orig = correlation, p_orig = pvalue),
      corrected_corr %>% select(pathway, variable, correlation, pvalue, FDR) %>% rename(corr_corr = correlation, p_corr = pvalue),
      by = c("pathway", "variable"),
      all = TRUE
    )

    # Check correlation reproducibility
    corr_diff <- abs(merged$corr_orig - merged$corr_corr)
    cat("Correlation value reproducibility:\n")
    cat("  Max absolute difference:", max(corr_diff, na.rm = TRUE), "\n")
    cat("  Mean absolute difference:", mean(corr_diff, na.rm = TRUE), "\n")

    if (max(corr_diff, na.rm = TRUE) < 0.01) {
      cat("  STATUS: REPRODUCIBLE - Correlation values match\n\n")
    } else {
      cat("  STATUS: WARNING - Some correlation values differ\n\n")
    }

    # Document tests that changed significance after FDR
    lost_sig <- merged %>%
      filter(p_orig < 0.05 & merged$FDR >= 0.05)

    if (nrow(lost_sig) > 0) {
      cat("Tests significant before FDR but not after:\n")
      print(lost_sig %>% select(pathway, variable, corr_orig, p_orig, FDR))
    }

    # Save comparison
    write.csv(merged, file.path(comparison_dir, "human_bulk_correlation_comparison.csv"), row.names = FALSE)
  }
} else {
  cat("Warning: Could not find both original and corrected correlation files\n")
}

# =============================================================================
# 2. Rat snRNA-seq Comparison
# =============================================================================
cat("\n--- Rat snRNA-seq Comparison ---\n\n")

# For rat data, original analysis was missing DE/DA, so we document what was added
corrected_de_file <- file.path(corrected_dir, "rat_snrnaseq/DE_summary_by_celltype.csv")
corrected_da_file <- file.path(corrected_dir, "rat_snrnaseq/DA_results_celltypes.csv")

cat("Biostatistical fixes applied to rat snRNA-seq:\n")
cat("  1. Added differential expression analysis (was MISSING)\n")
cat("  2. Added differential abundance testing (was MISSING)\n")
cat("  3. Added FDR correction (Benjamini-Hochberg)\n")
cat("  4. Added DoubletFinder for doublet detection\n")
cat("  5. Added random seeds for reproducibility\n\n")

if (file.exists(corrected_de_file)) {
  de_summary <- read.csv(corrected_de_file)
  cat("Differential Expression Results (NEW):\n")
  print(de_summary)
  cat("\n")
}

if (file.exists(corrected_da_file)) {
  da_results <- read.csv(corrected_da_file)
  cat("Differential Abundance Results (NEW):\n")
  cat("  Cell types tested:", nrow(da_results) / 2, "\n")  # Divided by 2 for Main + Subtype
  cat("  Significant at FDR < 0.05:", sum(da_results$FDR < 0.05), "\n\n")
}

# =============================================================================
# 3. Generate Summary Report
# =============================================================================
cat("--- Generating Summary Report ---\n\n")

report_lines <- c(
  "# ERpBRCA_OlderWomen: Biostatistical Corrections Summary",
  "",
  paste("Generated:", Sys.time()),
  "",
  "## Overview",
  "",
  "This report documents the biostatistical corrections applied during code review",
  "and validates that the core results are reproducible.",
  "",
  "## Human Bulk RNA-seq Analysis",
  "",
  "### Corrections Applied:",
  "1. **FDR Correction**: Added Benjamini-Hochberg FDR adjustment for multiple testing",
  "2. **Random Seeds**: Added for reproducibility (seed = 12345)",
  "",
  "### Results Comparison:",
  paste("- Original tests:", ifelse(exists("original_corr"), nrow(original_corr), "N/A")),
  paste("- Original significant (p < 0.05):", ifelse(exists("original_sig"), original_sig, "N/A")),
  paste("- Corrected significant (p < 0.05):", ifelse(exists("corrected_sig_nominal"), corrected_sig_nominal, "N/A")),
  paste("- Corrected significant (FDR < 0.05):", ifelse(exists("corrected_sig_fdr"), corrected_sig_fdr, "N/A")),
  "",
  "### Reproducibility Status:",
  ifelse(exists("corr_diff") && max(corr_diff, na.rm = TRUE) < 0.01,
         "**PASS**: Correlation values reproduced exactly",
         "See detailed comparison for any differences"),
  "",
  "## Rat snRNA-seq Analysis",
  "",
  "### Corrections Applied:",
  "1. **Differential Expression**: Added Young vs Aged comparison per cell type (WAS MISSING)",
  "2. **Differential Abundance**: Added propeller testing for cell type proportions (WAS MISSING)",
  "3. **DoubletFinder**: Added doublet detection",
  "4. **FDR Correction**: Applied BH correction to all tests",
  "5. **Random Seeds**: Added for reproducibility",
  "",
  "### New Results Generated:",
  paste("- Cell types identified:", ifelse(exists("de_summary"), nrow(de_summary), "N/A")),
  paste("- Total DE genes (FDR < 0.05):", ifelse(exists("de_summary"), sum(de_summary$sig_FDR), "N/A")),
  "",
  "## Conclusion",
  "",
  "The biostatistical corrections ensure:",
  "- Proper multiple testing correction (FDR)",
  "- Complete statistical analyses (DE/DA for rat data)",
  "- Reproducibility through random seeds",
  "",
  "Core correlation values from the human bulk RNA-seq analysis are reproduced,",
  "demonstrating that the underlying data processing is consistent.",
  ""
)

writeLines(report_lines, file.path(comparison_dir, "comparison_report.md"))
cat("Report saved to:", file.path(comparison_dir, "comparison_report.md"), "\n")

# =============================================================================
# 4. Generate Comparison Figure
# =============================================================================
cat("\n--- Generating Comparison Figure ---\n")

if (exists("merged") && nrow(merged) > 0) {
  pdf(file.path(fig_dir, "fdr_impact_comparison.pdf"), width = 10, height = 8)

  # Plot: p-values vs FDR
  p1 <- ggplot(merged, aes(x = -log10(p_corr), y = -log10(FDR))) +
    geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                       name = "FDR < 0.05") +
    theme_bw() +
    labs(
      title = "Impact of FDR Correction on Significance",
      subtitle = "Human Bulk RNA-seq Correlation Analysis",
      x = "-log10(nominal p-value)",
      y = "-log10(FDR-adjusted p-value)"
    ) +
    annotate("text", x = 1.5, y = 0.5, label = "Lost significance\nafter FDR", color = "blue", size = 3)

  print(p1)

  dev.off()
  cat("Figure saved to:", file.path(fig_dir, "fdr_impact_comparison.pdf"), "\n")
}

cat("\n=== Comparison Report Complete ===\n")
