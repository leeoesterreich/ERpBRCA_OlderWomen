#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/05_visualize.R
# Generate bubble plots and other visualizations
#
# Inputs:
#   - analysis/01_human_bulk_rnaseq/outputs/correlation_results.rds
#
# Outputs:
#   - figures/by_analysis/human_bulk_rnaseq/correlation_bubbleplot_original.pdf
#   - figures/by_analysis/human_bulk_rnaseq/correlation_bubbleplot_fdr.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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
fig_dir <- file.path(project_root, "figures/by_analysis/human_bulk_rnaseq")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Generating Visualizations ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Correlation Results
# -----------------------------------------------------------------------------
cat("Step 1: Loading correlation results...\n")
corr_file <- file.path(output_dir, "correlation_results.rds")
check_file_exists(corr_file, "Correlation results")
corr_results <- readRDS(corr_file)

# Filter out RAB19 and prepare for plotting
# Note: Document rationale - RAB19 was removed in original analysis
plot_data <- corr_results %>%
  filter(!grepl("RAB19", GeneSymb)) %>%
  mutate(
    GeneSymb = factor(GeneSymb, levels = c(
      "TFF1", "SAA1", "PGR", "PAK4", "HSD17B7", "HSD17B2",
      "GREB1", "ESR1", "CYP19A1", "Age"
    )),
    # Truncate long pathway names for display
    PathwayName_short = gsub("HALLMARK_", "HM_", PathwayName),
    PathwayName_short = gsub("REACTOME_", "RC_", PathwayName_short),
    PathwayName_short = gsub("GOBP_", "GO_", PathwayName_short)
  )

# Color palette (Wong colorblind-safe diverging: blue → white → vermillion)
my_palette <- colorRampPalette(c("#0072B2", "#F7F7F7", "#D55E00"))(100)

# -----------------------------------------------------------------------------
# Step 2: Original Bubble Plot (no FDR)
# -----------------------------------------------------------------------------
cat("Step 2: Creating original bubble plot...\n")

p_original <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName_short)) +
  geom_point(aes(size = -log10(Spearman_pval), color = Spearman_Rho)) +
  scale_color_gradientn("Spearman Rho", colors = my_palette, limits = c(-1, 1)) +
  scale_size_continuous("-log10(p)", range = c(1, 10)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.title = element_blank(),
    panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white")
  ) +
  coord_flip() +
  ggtitle("Gene-Pathway Correlations (Original, no FDR)")

# Save in multiple formats
ggsave(file.path(fig_dir, "correlation_bubbleplot_original.pdf"), p_original, width = 12, height = 7)
ggsave(file.path(fig_dir, "correlation_bubbleplot_original.png"), p_original, width = 12, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# Step 3: FDR-Corrected Bubble Plot
# -----------------------------------------------------------------------------
cat("Step 3: Creating FDR-corrected bubble plot...\n")

p_fdr <- ggplot(plot_data, aes(x = GeneSymb, y = PathwayName_short)) +
  geom_point(aes(size = -log10(FDR_qval), color = Spearman_Rho)) +
  scale_color_gradientn("Spearman Rho", colors = my_palette, limits = c(-1, 1)) +
  scale_size_continuous("-log10(FDR)", range = c(1, 10)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.title = element_blank(),
    panel.border = element_rect(linewidth = 0.7, linetype = "solid", colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white")
  ) +
  coord_flip() +
  ggtitle("Gene-Pathway Correlations (FDR-Corrected)")

# Save in multiple formats
ggsave(file.path(fig_dir, "correlation_bubbleplot_fdr.pdf"), p_fdr, width = 12, height = 7)
ggsave(file.path(fig_dir, "correlation_bubbleplot_fdr.png"), p_fdr, width = 12, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# Step 4: Side-by-side Comparison
# -----------------------------------------------------------------------------
cat("Step 4: Creating side-by-side comparison...\n")

p_combined <- p_original + p_fdr +
  plot_annotation(
    title = "Effect of FDR Correction on Gene-Pathway Correlations",
    subtitle = paste0(
      "Significant at p<0.05: ", sum(plot_data$Sig_nominal), " | ",
      "Significant at FDR<0.05: ", sum(plot_data$Sig_FDR)
    )
  )

# Save in multiple formats
ggsave(file.path(fig_dir, "correlation_bubbleplot_comparison.pdf"), p_combined, width = 20, height = 8)
ggsave(file.path(fig_dir, "correlation_bubbleplot_comparison.png"), p_combined, width = 20, height = 8, dpi = 300)

cat("\n=== Visualization complete ===\n")
cat("Figures saved to:", fig_dir, "\n")
