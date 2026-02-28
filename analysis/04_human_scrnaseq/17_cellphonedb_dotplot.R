#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/17_cellphonedb_dotplot.R
# Generate Figure 7D - CellPhoneDB comparative dot plot (Young vs Elderly)
#
# Inputs:
#   - outputs/cellphonedb/results_elderly/pvalues.csv
#   - outputs/cellphonedb/results_elderly/means.csv
#   - outputs/cellphonedb/results_young/pvalues.csv
#   - outputs/cellphonedb/results_young/means.csv
#
# Outputs:
#   - figures/fig7d_cellphonedb_dotplot.png

set.seed(12345)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(RColorBrewer)
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
cpdb_dir <- file.path(script_dir, "outputs/cellphonedb")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== CellPhoneDB Comparative Dot Plot (Figure 7D) ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load CellPhoneDB results for both age groups
# -----------------------------------------------------------------------------
cat("Step 1: Loading CellPhoneDB results...\n")

load_cpdb_results <- function(group) {
  results_dir <- file.path(cpdb_dir, paste0("results_", tolower(group)))

  # CellPhoneDB v5 uses timestamped file names
  pvals_files <- list.files(results_dir, pattern = "pvalues.*\\.txt$", full.names = TRUE)
  means_files <- list.files(results_dir, pattern = "means.*\\.txt$", full.names = TRUE)
  # Filter out significant_means
  means_files <- means_files[!grepl("significant", means_files)]

  if (length(pvals_files) == 0) {
    pvals_file <- file.path(results_dir, "pvalues.csv")
    means_file <- file.path(results_dir, "means.csv")
    if (!file.exists(pvals_file)) {
      cat(sprintf("  WARNING: No pvalues file found in %s\n", results_dir))
      return(NULL)
    }
  } else {
    pvals_file <- sort(pvals_files, decreasing = TRUE)[1]
    means_file <- sort(means_files, decreasing = TRUE)[1]
  }

  cat(sprintf("  Loading %s: %s\n", group, basename(pvals_file)))

  pvals <- fread(pvals_file)
  means <- fread(means_file)

  list(pvals = pvals, means = means, group = group)
}

results_elderly <- load_cpdb_results("Elderly")
results_young <- load_cpdb_results("Young")

if (is.null(results_elderly) || is.null(results_young)) {
  stop("Need both Young and Elderly CellPhoneDB results for comparative analysis.")
}

# -----------------------------------------------------------------------------
# Step 2: Process both datasets
# -----------------------------------------------------------------------------
cat("\nStep 2: Processing both age groups...\n")

process_cpdb <- function(results) {
  pvals <- results$pvals
  means <- results$means
  group <- results$group

  # Find cell pair columns (columns with "|" in name)
  cell_pair_cols <- colnames(pvals)[grepl("\\|", colnames(pvals))]

  # Convert to long format
  pvals_long <- pvals %>%
    select(interacting_pair, all_of(cell_pair_cols)) %>%
    pivot_longer(cols = -interacting_pair,
                 names_to = "cell_pair",
                 values_to = "pvalue")

  means_long <- means %>%
    select(interacting_pair, all_of(cell_pair_cols)) %>%
    pivot_longer(cols = -interacting_pair,
                 names_to = "cell_pair",
                 values_to = "mean_expr")

  # Merge and add group
  combined <- left_join(pvals_long, means_long,
                        by = c("interacting_pair", "cell_pair")) %>%
    mutate(age_group = group)

  combined
}

data_elderly <- process_cpdb(results_elderly)
data_young <- process_cpdb(results_young)

# Combine both groups
all_data <- bind_rows(data_young, data_elderly)
cat(sprintf("  Combined data: %d rows\n", nrow(all_data)))

# -----------------------------------------------------------------------------
# Step 3: Identify significant interactions in either group
# -----------------------------------------------------------------------------
cat("\nStep 3: Identifying significant interactions...\n")

# Focus on MACROPHAGE communication with ADAPTIVE IMMUNE cells
# Per manuscript legend: "macrophages from older patients have less communication
# with adaptive immune cells compared to macrophages from younger patients"
adaptive_immune <- c("T cells", "Tcells", "Bcells", "NK", "NKT", "Plasma")

# Filter to Macrophage | AdaptiveImmune OR AdaptiveImmune | Macrophage pairs
immune_data <- all_data %>%
  filter(
    (grepl("Macrophage", cell_pair, ignore.case = TRUE) &
     grepl(paste(adaptive_immune, collapse = "|"), cell_pair, ignore.case = TRUE))
  )

cat(sprintf("  Immune-relevant: %d rows\n", nrow(immune_data)))

# Find L-R pairs significant in at least one age group
sig_lr_pairs <- immune_data %>%
  filter(!is.na(pvalue) & pvalue < 0.05) %>%
  group_by(interacting_pair) %>%
  summarise(
    n_sig = n(),
    n_young_sig = sum(age_group == "Young"),
    n_elderly_sig = sum(age_group == "Elderly"),
    min_pval = min(pvalue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # For macrophage-focused analysis, include pairs significant in either group
  filter(n_sig >= 1) %>%
  arrange(min_pval) %>%
  head(50)  # Top 50 L-R pairs

cat(sprintf("  Significant L-R pairs: %d\n", nrow(sig_lr_pairs)))

# Find top cell pairs with most interactions
sig_cell_pairs <- immune_data %>%
  filter(!is.na(pvalue) & pvalue < 0.05,
         interacting_pair %in% sig_lr_pairs$interacting_pair) %>%
  group_by(cell_pair) %>%
  summarise(n_interactions = n(), .groups = "drop") %>%
  arrange(desc(n_interactions)) %>%
  head(15)  # Top 15 cell pairs

cat(sprintf("  Top cell pairs: %d\n", nrow(sig_cell_pairs)))

# Create final plot data
plot_data <- immune_data %>%
  filter(interacting_pair %in% sig_lr_pairs$interacting_pair,
         cell_pair %in% sig_cell_pairs$cell_pair) %>%
  mutate(
    neg_log_pval = -log10(pvalue + 1e-10),
    neg_log_pval = pmin(neg_log_pval, 3),  # Cap at 3 to match manuscript scale
    # Use log2 mean expression to match manuscript
    log2_mean = log2(mean_expr + 1),
    # Set non-significant to NA for cleaner plot
    log2_mean_plot = ifelse(pvalue < 0.05, log2_mean, NA)
  )

cat(sprintf("  Final plot data: %d points\n", nrow(plot_data)))

# -----------------------------------------------------------------------------
# Step 4: Generate comparative dot plot
# -----------------------------------------------------------------------------
cat("\nStep 4: Generating comparative dot plot...\n")

# Rename age groups to match manuscript convention and embed in cell_pair
plot_data <- plot_data %>%
  mutate(
    age_label = case_when(
      age_group == "Young" ~ "Younger",
      age_group == "Elderly" ~ "Older",
      TRUE ~ age_group
    ),
    # Create cell_pair_age: simpler format "cell_pair (Age)"
    cell_pair_age = paste0(cell_pair, " (", age_label, ")")
  )

# Order L-R pairs by mean expression (simpler approach to avoid pivot issues)
lr_order <- plot_data %>%
  filter(!is.na(log2_mean_plot)) %>%
  group_by(interacting_pair) %>%
  summarise(mean_log2 = mean(log2_mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_log2) %>%
  pull(interacting_pair)

plot_data$interacting_pair <- factor(plot_data$interacting_pair, levels = lr_order)

# Order cell pairs: group by base pair, then by age (Younger before Older)
plot_data <- plot_data %>%
  arrange(cell_pair, age_label)

# Create ordered factor for cell_pair_age
cell_pair_age_order <- unique(plot_data$cell_pair_age)
plot_data$cell_pair_age <- factor(plot_data$cell_pair_age, levels = cell_pair_age_order)

# Create single-panel dot plot (no faceting) - matches manuscript structure
p <- ggplot(plot_data, aes(x = cell_pair_age, y = interacting_pair)) +
  geom_point(aes(size = neg_log_pval, color = log2_mean_plot), na.rm = TRUE) +
  scale_size_continuous(
    name = "-log10(p-value)",
    range = c(0.5, 5),
    breaks = c(1, 2, 3),
    limits = c(0, 3)
  ) +
  scale_color_gradientn(
    name = "Log2 mean",
    colors = c("blue", "yellow", "red"),
    na.value = "grey90",
    limits = c(-1, 5)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  labs(
    x = "Cell type pairs",
    y = "Ligand-Receptor pairs",
    title = "CellPhoneDB: Macrophage-Adaptive Immune Communication"
  )

# Save
output_file <- file.path(figures_dir, "fig7d_cellphonedb_dotplot.png")
ggsave(output_file, p, width = 14, height = 12, dpi = 300)
cat(sprintf("  Saved: %s\n", output_file))

# Also save as PDF
pdf_file <- file.path(figures_dir, "fig7d_cellphonedb_dotplot.pdf")
ggsave(pdf_file, p, width = 14, height = 12)
cat(sprintf("  Saved: %s\n", pdf_file))

# -----------------------------------------------------------------------------
# Step 5: Save summary statistics
# -----------------------------------------------------------------------------
cat("\nStep 5: Saving summary statistics...\n")

summary_df <- plot_data %>%
  filter(!is.na(log2_mean_plot)) %>%
  group_by(interacting_pair, age_label) %>%
  summarise(
    n_cell_pairs = n(),
    mean_pval = mean(pvalue, na.rm = TRUE),
    min_pval = min(pvalue, na.rm = TRUE),
    log2_mean = mean(log2_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = age_label,
    values_from = c(n_cell_pairs, mean_pval, min_pval, log2_mean),
    names_sep = "_"
  ) %>%
  arrange(pmin(min_pval_Younger, min_pval_Older, na.rm = TRUE))

output_csv <- file.path(script_dir, "outputs", "cellphonedb_summary.csv")
fwrite(summary_df, output_csv)
cat(sprintf("  Saved: %s\n", output_csv))

cat("\n=== CellPhoneDB comparative dot plot complete ===\n")
