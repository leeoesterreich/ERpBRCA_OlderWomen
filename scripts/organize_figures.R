#!/usr/bin/env Rscript
# scripts/organize_figures.R
# Copy figures from analysis outputs to manuscript folder with consistent naming
#
# Outputs:
#   - figures/manuscript/Figure_1_*.pdf
#   - figures/manuscript/Figure_2_*.pdf
#   - etc.

set.seed(12345)

script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, ".."))

# Create manuscript figures directory
manuscript_dir <- file.path(project_root, "figures/manuscript")
dir.create(manuscript_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Organizing Figures for Manuscript ===\n")

# Define figure mappings: source -> destination
# These mappings should be updated based on manuscript requirements
figure_map <- list(
  # Figure 1: Human bulk RNA-seq overview
  list(
    src = "figures/by_analysis/01_human_bulk_rnaseq/pca_plot.pdf",
    dst = "Figure_1a_PCA_samples.pdf"
  ),
  list(
    src = "figures/by_analysis/01_human_bulk_rnaseq/gsva_heatmap.pdf",
    dst = "Figure_1b_GSVA_heatmap.pdf"
  ),
  list(
    src = "figures/by_analysis/01_human_bulk_rnaseq/correlation_bubbleplot.pdf",
    dst = "Figure_1c_correlation_bubbleplot.pdf"
  ),

  # Figure 2: Rat snRNA-seq overview
  list(
    src = "figures/by_analysis/02_rat_snrnaseq/umap_celltypes.pdf",
    dst = "Figure_2a_UMAP_celltypes.pdf"
  ),
  list(
    src = "figures/by_analysis/02_rat_snrnaseq/cell_proportions.pdf",
    dst = "Figure_2b_cell_proportions.pdf"
  ),

  # Figure 3: Immune analysis
  list(
    src = "figures/by_analysis/02_rat_snrnaseq/myeloid_umap.pdf",
    dst = "Figure_3a_myeloid_UMAP.pdf"
  ),
  list(
    src = "figures/by_analysis/02_rat_snrnaseq/nkt_umap.pdf",
    dst = "Figure_3b_NKT_UMAP.pdf"
  ),

  # Supplementary figures
  list(
    src = "figures/by_analysis/01_human_bulk_rnaseq/qc_library_sizes.pdf",
    dst = "Supp_Figure_1_library_sizes.pdf"
  ),
  list(
    src = "figures/by_analysis/02_rat_snrnaseq/qc_violin_plots.pdf",
    dst = "Supp_Figure_2_QC_violins.pdf"
  ),

  # Comparison figures
  list(
    src = "results/comparison/figures/bubbleplot_sidebyside.pdf",
    dst = "Supp_Figure_comparison_bubbleplots.pdf"
  ),
  list(
    src = "results/comparison/figures/significance_change_heatmap.pdf",
    dst = "Supp_Figure_significance_changes.pdf"
  )
)

# Copy files
n_copied <- 0
n_missing <- 0

for (mapping in figure_map) {
  src_path <- file.path(project_root, mapping$src)
  dst_path <- file.path(manuscript_dir, mapping$dst)

  if (file.exists(src_path)) {
    file.copy(src_path, dst_path, overwrite = TRUE)
    cat("  Copied:", mapping$dst, "\n")
    n_copied <- n_copied + 1
  } else {
    cat("  Missing:", mapping$src, "\n")
    n_missing <- n_missing + 1
  }
}

cat("\n=== Figure organization complete ===\n")
cat("Copied:", n_copied, "figures\n")
cat("Missing:", n_missing, "source files\n")
cat("Output directory:", manuscript_dir, "\n")
