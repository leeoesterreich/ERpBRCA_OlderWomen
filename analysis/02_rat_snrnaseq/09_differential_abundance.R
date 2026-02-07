#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/09_differential_abundance.R
# Differential abundance analysis of cell types between Young vs Aged
# BIOSTATISTICAL FIX: This statistical test was missing from original
#
# Inputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - results/corrected/rat_snrnaseq/DA_results_celltypes.csv
#   - figures/by_analysis/rat_snrnaseq/DA_proportion_plots.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(speckle)  # For propeller differential abundance
  library(tidyr)
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
results_dir <- file.path(project_root, "results/corrected/rat_snrnaseq")
fig_dir <- file.path(project_root, "figures/by_analysis/rat_snrnaseq")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Differential Abundance Analysis ===\n")
cat("Random seed: 12345\n")
cat("Project root:", project_root, "\n")
cat("Results directory:", results_dir, "\n\n")

# Input file
input_rds <- file.path(output_dir, "seurat_annotated.rds")

# Verify input file exists
cat("Checking input files...\n")
check_file_exists(input_rds, "Annotated Seurat object")
cat("\n")

# -----------------------------------------------------------------------------
# Step 1: Load Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")
seurat_obj <- readRDS(input_rds)

# Validate Seurat object
if (!inherits(seurat_obj, "Seurat")) {
  stop("ERROR: Loaded object is not a valid Seurat object")
}
if (ncol(seurat_obj) == 0) {
  stop("ERROR: Seurat object contains no cells")
}

# Get metadata
meta <- seurat_obj@meta.data %>%
  select(orig.ident, AgeGroup, CellTypeByMarker_RatsnRNAseq, CellTypeMacroTcell_RatsnRNAseq)

cat("  Total cells:", nrow(meta), "\n")
cat("  Samples:", length(unique(meta$orig.ident)), "\n")
cat("  Age groups:\n")
print(table(meta$AgeGroup))

# Check sample sizes - document limitation
n_samples <- meta %>%
  select(orig.ident, AgeGroup) %>%
  distinct()

cat("\n  Sample distribution:\n")
print(table(n_samples$AgeGroup))

# -----------------------------------------------------------------------------
# Step 2: Calculate Cell Type Proportions
# -----------------------------------------------------------------------------
cat("\nStep 2: Calculating cell type proportions...\n")

# Main cell types
prop_main <- meta %>%
  group_by(orig.ident, AgeGroup, CellTypeByMarker_RatsnRNAseq) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  )

cat("  Cell type proportions calculated\n")
cat("  Cell types:", length(unique(prop_main$CellTypeByMarker_RatsnRNAseq)), "\n")

# Detailed subtypes (if different from main)
prop_sub <- meta %>%
  group_by(orig.ident, AgeGroup, CellTypeMacroTcell_RatsnRNAseq) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  )

# -----------------------------------------------------------------------------
# Step 3: Run Propeller DA Test
# -----------------------------------------------------------------------------
cat("\nStep 3: Running propeller differential abundance test...\n")

# Prepare data for propeller
# Note: propeller requires counts, sample IDs, and cluster IDs

# Check for valid cell types before running propeller
cell_types <- unique(seurat_obj$CellTypeByMarker_RatsnRNAseq)
cell_types <- cell_types[!is.na(cell_types)]
cat("  Cell types for DA:", length(cell_types), "-", paste(cell_types, collapse = ", "), "\n")

if (length(cell_types) < 2) {
  warning("Need at least 2 cell types for differential abundance analysis")
  cat("\nWARNING: Fewer than 2 cell types found. Skipping propeller analysis.\n")
  cat("This may indicate cell type annotation issues.\n")

  # Still generate proportion plots with available data
  da_main <- NULL
} else {
  # Main cell types
  da_main <- NULL
  tryCatch({
    da_main <- propeller(
      clusters = seurat_obj$CellTypeByMarker_RatsnRNAseq,
      sample = seurat_obj$orig.ident,
      group = seurat_obj$AgeGroup
    )

    da_main$CellType <- rownames(da_main)
    # BIOSTATISTICAL FIX: Apply BH FDR correction
    da_main$FDR <- p.adjust(da_main$P.Value, method = "BH")
    da_main$Level <- "Main"

    cat("\nMain cell type results:\n")
    print(da_main %>% select(CellType, PropMean.Aged, PropMean.Young, P.Value, FDR) %>% arrange(P.Value))

  }, error = function(e) {
    cat("  Error in main cell type propeller:", e$message, "\n")
  })
}

# Detailed subtypes (if enough cells per category)
da_sub <- NULL
tryCatch({
  da_sub <- propeller(
    clusters = seurat_obj$CellTypeMacroTcell_RatsnRNAseq,
    sample = seurat_obj$orig.ident,
    group = seurat_obj$AgeGroup
  )

  da_sub$CellType <- rownames(da_sub)
  da_sub$FDR <- p.adjust(da_sub$P.Value, method = "BH")
  da_sub$Level <- "Subtype"

  cat("\nSubtype results:\n")
  print(da_sub %>% select(CellType, PropMean.Aged, PropMean.Young, P.Value, FDR) %>% arrange(P.Value))

}, error = function(e) {
  cat("  Warning: Subtype analysis failed -", e$message, "\n")
})

# Combine results
da_all <- NULL
if (!is.null(da_main) && !is.null(da_sub)) {
  da_all <- rbind(da_main, da_sub)
} else if (!is.null(da_main)) {
  da_all <- da_main
} else if (!is.null(da_sub)) {
  da_all <- da_sub
}

# -----------------------------------------------------------------------------
# Step 4: Generate Proportion Plots
# -----------------------------------------------------------------------------
cat("\nStep 4: Generating proportion plots...\n")

pdf(file.path(fig_dir, "DA_proportion_plots.pdf"), width = 14, height = 10)

# Stacked bar plot by sample
p1 <- ggplot(prop_main, aes(x = orig.ident, y = proportion, fill = CellTypeByMarker_RatsnRNAseq)) +
  geom_bar(stat = "identity") +
  facet_wrap(~AgeGroup, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Cell Type Proportions by Sample",
    x = "Sample", y = "Proportion", fill = "Cell Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)

# Box plot comparing proportions
p2 <- ggplot(prop_main, aes(x = CellTypeByMarker_RatsnRNAseq, y = proportion * 100, fill = AgeGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.7, size = 2) +
  theme_bw() +
  labs(
    title = "Cell Type Proportions: Young vs Aged",
    subtitle = "Note: n=3 per group - interpret with caution",
    x = "Cell Type", y = "Percentage", fill = "Age Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p2)

# DA results visualization (only if we have results)
if (!is.null(da_all)) {
  da_plot_data <- da_all %>%
    mutate(
      Sig = FDR < 0.05,
      LogFC = log2(PropMean.Aged / PropMean.Young)
    )

  # Handle infinite values (log how many are affected)
  n_infinite <- sum(is.infinite(da_plot_data$LogFC))
  if (n_infinite > 0) {
    cat("  Warning:", n_infinite, "infinite LogFC values set to NA (division by zero)\n")
  }
  da_plot_data$LogFC[is.infinite(da_plot_data$LogFC)] <- NA

  p3 <- ggplot(da_plot_data %>% filter(!is.na(LogFC)),
               aes(x = reorder(CellType, LogFC), y = LogFC, fill = Sig)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    facet_wrap(~Level, scales = "free_y") +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    theme_bw() +
    labs(
      title = "Differential Abundance: Aged vs Young",
      subtitle = "Red = FDR < 0.05",
      x = "Cell Type", y = "log2 Fold Change (Aged/Young)"
    )
  print(p3)

  # P-value bar plot
  p4 <- ggplot(da_plot_data, aes(x = reorder(CellType, -log10(FDR)), y = -log10(FDR), fill = Sig)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    coord_flip() +
    facet_wrap(~Level, scales = "free_y") +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    theme_bw() +
    labs(
      title = "Differential Abundance Significance",
      subtitle = "Dashed line = FDR 0.05",
      x = "Cell Type", y = "-log10(FDR)"
    )
  print(p4)
} else {
  cat("  Note: Skipping DA plots - no propeller results available\n")
}

dev.off()

cat("  Proportion plots saved to:", file.path(fig_dir, "DA_proportion_plots.pdf"), "\n")

# -----------------------------------------------------------------------------
# Step 5: Save Results
# -----------------------------------------------------------------------------
cat("\nStep 5: Saving results...\n")

# Save proportions regardless of DA results
write.csv(prop_main, file.path(results_dir, "cell_proportions_by_sample.csv"), row.names = FALSE)

# Save DA results if available
if (!is.null(da_all)) {
  # Add sample size warning to results
  da_all$Warning <- "n=3 per group; interpret with caution"
  write.csv(da_all, file.path(results_dir, "DA_results_celltypes.csv"), row.names = FALSE)
} else {
  cat("  Note: No DA results to save (fewer than 2 cell types)\n")
}

# Also document the limitation
writeLines(
  c(
    "# Differential Abundance Analysis Notes",
    "",
    "## Sample Size Limitation",
    "- Young: n=3 samples",
    "- Aged: n=3 samples",
    "",
    "With only 3 samples per group, statistical power is severely limited.",
    "Effect sizes and biological trends may be more meaningful than p-values.",
    "",
    "## Method",
    "- propeller (speckle package) for differential abundance",
    "- Benjamini-Hochberg FDR correction",
    "",
    "## Interpretation Guidelines",
    "1. Focus on effect sizes (log2 fold changes) rather than p-values alone",
    "2. Consider biological plausibility of observed changes",
    "3. Validate findings with independent data or methods if possible",
    "4. Report results with appropriate caveats about sample size",
    ""
  ),
  file.path(results_dir, "DA_analysis_notes.md")
)

cat("  Results saved:\n")
cat("    - DA_results_celltypes.csv\n")
cat("    - cell_proportions_by_sample.csv\n")
cat("    - DA_analysis_notes.md (limitations documentation)\n")

cat("\n=== Differential abundance analysis complete ===\n")
cat("\nIMPORTANT: With n=3 per group, statistical power is limited.\n")
cat("Interpret p-values with caution; effect sizes may be more informative.\n")
