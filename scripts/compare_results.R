#!/usr/bin/env Rscript
# scripts/compare_results.R
# Compare original vs corrected results and generate comparison report
#
# Outputs:
#   - results/comparison/correlation_comparison.csv
#   - results/comparison/significant_findings_summary.md
#   - results/comparison/comparison_report.html
#   - results/comparison/comparison_report.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(rmarkdown)
  library(knitr)
})

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, ".."))
results_dir <- file.path(project_root, "results/comparison")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Generating Comparison Report ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Correlation Results
# -----------------------------------------------------------------------------
cat("Step 1: Loading correlation results...\n")

original_file <- file.path(project_root, "results/original/human_bulk_rnaseq/correlations_no_fdr.csv")
corrected_file <- file.path(project_root, "results/corrected/human_bulk_rnaseq/correlations_with_fdr.csv")

if (!file.exists(original_file) || !file.exists(corrected_file)) {
  cat("Warning: Correlation files not found. Run analyses first.\n")
  original <- NULL
  corrected <- NULL
} else {
  original <- read.csv(original_file)
  corrected <- read.csv(corrected_file)
}

# -----------------------------------------------------------------------------
# Step 2: Create Comparison Table
# -----------------------------------------------------------------------------
cat("Step 2: Creating comparison table...\n")

if (!is.null(corrected)) {
  comparison <- corrected %>%
    mutate(
      Status = case_when(
        Sig_nominal & Sig_FDR ~ "Remains Significant",
        Sig_nominal & !Sig_FDR ~ "Lost Significance (FDR)",
        !Sig_nominal ~ "Not Significant"
      )
    ) %>%
    select(GeneSymb, PathwayName, Spearman_Rho, Spearman_pval, FDR_qval, Status)

  # Summary statistics
  n_total <- nrow(comparison)
  n_sig_original <- sum(comparison$Spearman_pval < 0.05)
  n_sig_fdr <- sum(comparison$FDR_qval < 0.05)
  n_lost <- sum(comparison$Status == "Lost Significance (FDR)")

  cat("  Total tests:", n_total, "\n")
  cat("  Significant (p<0.05):", n_sig_original, "\n")
  cat("  Significant (FDR<0.05):", n_sig_fdr, "\n")
  cat("  Lost significance:", n_lost, "\n")

  write.csv(comparison, file.path(results_dir, "correlation_comparison.csv"), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# Step 3: Generate Summary Markdown
# -----------------------------------------------------------------------------
cat("Step 3: Generating summary markdown...\n")

summary_md <- c(
  "# Biostatistical Corrections: Results Comparison",
  "",
  "## Overview",
  "",
  "This document summarizes the impact of biostatistical corrections on the analysis results.",
  "",
  "## Human Bulk RNA-seq Correlations",
  "",
  if (!is.null(corrected)) {
    c(
      paste("- **Total correlation tests:**", n_total),
      paste("- **Significant at p < 0.05 (original):**", n_sig_original),
      paste("- **Significant at FDR < 0.05 (corrected):**", n_sig_fdr),
      paste("- **Lost significance after FDR:**", n_lost),
      "",
      "### Key Findings Status",
      ""
    )
  } else {
    "Results not yet available. Run the bulk RNA-seq analysis first."
  },
  "",
  "## Rat snRNA-seq Analysis",
  "",
  "### New Analyses Added",
  "- Differential expression (Young vs Aged) per cell type",
  "- Differential abundance testing using propeller",
  "- Doublet detection enabled",
  "",
  "### Sample Size Limitation",
  "- n=3 per age group",
  "- Statistical power is limited",
  "- Effect sizes may be more informative than p-values",
  "",
  "## Corrections Applied",
  "",
  "### Bulk RNA-seq",

  "1. DESeq2 design: `~1` -> `~ AgeRange + Group`",
  "2. Multiple testing: None -> Benjamini-Hochberg FDR",
  "3. Random seeds added for PROGENy",
  "4. QC plots added (PCA, library size, sample distances)",
  "",
  "### snRNA-seq",
  "1. DoubletFinder enabled (was commented out)",
  "2. Random seeds added for UMAP, Harmony, clustering",
  "3. Differential expression analysis added",
  "4. Differential abundance analysis added",
  "5. Cell type validation with reference datasets"
)

writeLines(summary_md, file.path(results_dir, "significant_findings_summary.md"))

# -----------------------------------------------------------------------------
# Step 4: Generate Comparison Figures
# -----------------------------------------------------------------------------
cat("Step 4: Generating comparison figures...\n")

fig_dir <- file.path(results_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

if (!is.null(corrected)) {
  # Bubble plot comparison
  pdf(file.path(fig_dir, "bubbleplot_sidebyside.pdf"), width = 16, height = 8)

  plot_data <- corrected %>%
    filter(!grepl("RAB19", GeneSymb)) %>%
    mutate(
      GeneSymb = factor(GeneSymb, levels = c(
        "TFF1", "SAA1", "PGR", "PAK4", "HSD17B7", "HSD17B2",
        "GREB1", "ESR1", "CYP19A1", "Age"
      ))
    )

  my_palette <- colorRampPalette(c("blue", "dodgerblue", "yellow", "orange", "red"))(100)

  p1 <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName)) +
    geom_point(aes(size = -log10(Spearman_pval), color = Spearman_Rho)) +
    scale_color_gradientn("Rho", colors = my_palette, limits = c(-1, 1)) +
    scale_size_continuous("-log10(p)", range = c(1, 8)) +
    coord_flip() + theme_bw() +
    ggtitle("Original (no FDR)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  p2 <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName)) +
    geom_point(aes(size = -log10(FDR_qval), color = Spearman_Rho)) +
    scale_color_gradientn("Rho", colors = my_palette, limits = c(-1, 1)) +
    scale_size_continuous("-log10(FDR)", range = c(1, 8)) +
    coord_flip() + theme_bw() +
    ggtitle("Corrected (with FDR)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  print(p1 + p2 + plot_annotation(title = "Effect of FDR Correction"))
  dev.off()

  # Significance change plot
  pdf(file.path(fig_dir, "significance_change_heatmap.pdf"), width = 10, height = 8)

  p3 <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName, fill = Status)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c(
      "Remains Significant" = "#2ECC71",
      "Lost Significance (FDR)" = "#E74C3C",
      "Not Significant" = "#BDC3C7"
    )) +
    coord_flip() + theme_bw() +
    ggtitle("Significance Status After FDR Correction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print(p3)
  dev.off()
}

# -----------------------------------------------------------------------------
# Step 5: Generate HTML and PDF Report
# -----------------------------------------------------------------------------
cat("Step 5: Generating HTML and PDF reports...\n")

# Create R Markdown content
rmd_content <- c(
  "---",
  'title: "Biostatistical Corrections: Comparison Report"',
  'author: "Automated Analysis"',
  paste0('date: "', Sys.Date(), '"'),
  "output:",
  "  html_document:",
  "    toc: true",
  "    toc_float: true",
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
  "```",
  "",
  readLines(file.path(results_dir, "significant_findings_summary.md")),
  "",
  "## Figures",
  "",
  "### Correlation Bubble Plots",
  "",
  paste0("![Side-by-side comparison](figures/bubbleplot_sidebyside.pdf){width=100%}"),
  "",
  "### Significance Changes",
  "",
  paste0("![Significance status](figures/significance_change_heatmap.pdf){width=80%}")
)

rmd_file <- file.path(results_dir, "comparison_report.Rmd")
writeLines(rmd_content, rmd_file)

# Render to HTML
tryCatch({
  rmarkdown::render(
    rmd_file,
    output_format = "html_document",
    output_file = "comparison_report.html",
    quiet = TRUE
  )
  cat("  HTML report generated\n")
}, error = function(e) {
  cat("  Warning: HTML generation failed -", e$message, "\n")
})

# Convert to PDF using pagedown
tryCatch({
  pagedown::chrome_print(
    file.path(results_dir, "comparison_report.html"),
    output = file.path(results_dir, "comparison_report.pdf")
  )
  cat("  PDF report generated\n")
}, error = function(e) {
  cat("  Warning: PDF generation failed -", e$message, "\n")
  cat("  You may need to install chromium: conda install -c conda-forge chromium\n")
})

cat("\n=== Comparison report complete ===\n")
cat("Outputs in:", results_dir, "\n")
