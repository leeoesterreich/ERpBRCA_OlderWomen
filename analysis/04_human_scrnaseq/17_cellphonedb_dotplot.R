#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/17_cellphonedb_dotplot.R
# Generate Figure 7D - CellPhoneDB comparative dot plot (Young vs Elderly)
#
# Usage:
#   Rscript 17_cellphonedb_dotplot.R [--mode=top50|curated]
#
# Modes:
#   top50   - Statistical filtering: top 50 L-R pairs by p-value (default)
#   curated - Use curated cytokine/chemokine list from original manuscript
#
# Inputs:
#   - outputs/cellphonedb/results_elderly/pvalues.csv
#   - outputs/cellphonedb/results_elderly/means.csv
#   - outputs/cellphonedb/results_young/pvalues.csv
#   - outputs/cellphonedb/results_young/means.csv
#   - outputs/curated_lr_pairs.txt (for curated mode)
#
# Outputs:
#   - figures/fig7d_cellphonedb_dotplot.png (or _curated.png)

set.seed(12345)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(RColorBrewer)
})

# -----------------------------------------------------------------------------
# Parse command line arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
filter_mode <- "top50"  # default

for (arg in args) {
  if (grepl("^--mode=", arg)) {
    filter_mode <- sub("^--mode=", "", arg)
  }
}

if (!filter_mode %in% c("top50", "curated")) {
  stop("Invalid mode. Use --mode=top50 or --mode=curated")
}

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
cat(sprintf("Filter mode: %s\n", filter_mode))

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

# Apply BH-FDR correction across all L-R pair tests
all_data <- all_data %>%
  mutate(pvalue_fdr = p.adjust(pvalue, method = "BH"))

cat(sprintf("  Combined data: %d rows\n", nrow(all_data)))
cat(sprintf("  FDR-corrected p-values applied (BH method)\n"))

# -----------------------------------------------------------------------------
# Step 3: Identify significant interactions in either group
# -----------------------------------------------------------------------------
cat("\nStep 3: Identifying significant interactions...\n")

# Focus on MACROPHAGE as SENDER only (manuscript Figure 7D)
# Exact target names from CellPhoneDB output (with spaces)
target_pairs <- c(
  "Macrophage|Bcells",
  "Macrophage|Cycling Tcells",
  "Macrophage|NK cells",
  "Macrophage|NKT cells",
  "Macrophage|T cells CD4",
  "Macrophage|T cells CD8"
)

# Filter to exact Macrophage|Target pairs
immune_data <- all_data %>%
  filter(cell_pair %in% target_pairs)

cat(sprintf("  Immune-relevant: %d rows\n", nrow(immune_data)))

# -----------------------------------------------------------------------------
# L-R pair selection based on filter_mode
# -----------------------------------------------------------------------------
if (filter_mode == "curated") {
  # Load curated L-R pairs (CellPhoneDB-compatible names)
  curated_file <- file.path(script_dir, "outputs/curated_lr_pairs.txt")
  if (!file.exists(curated_file)) {
    stop("Curated L-R pairs file not found: ", curated_file)
  }

  curated_pairs <- readLines(curated_file)
  # Remove comments and empty lines
  curated_pairs <- curated_pairs[!grepl("^#", curated_pairs) & nchar(trimws(curated_pairs)) > 0]
  curated_pairs <- trimws(curated_pairs)

  cat(sprintf("  Loaded %d curated L-R pairs\n", length(curated_pairs)))

  # FIX: Normalize L-R pair names to handle ordering differences (CD74_APP vs APP_CD74)
  # CellPhoneDB v5 may output pairs in different order than manuscript naming
  normalize_lr_pair <- function(pair) {
    # Handle complex names: remove common prefixes, normalize delimiters
    pair <- gsub("integrin_", "", pair, ignore.case = TRUE)
    pair <- gsub(" ", "_", pair)
    pair <- toupper(pair)
    # Sort components alphabetically to normalize A_B vs B_A
    parts <- strsplit(pair, "_")[[1]]
    if (length(parts) == 2) {
      return(paste(sort(parts), collapse = "_"))
    }
    # For complex pairs (3+ parts), keep as-is
    return(pair)
  }

  # Create normalized lookup table for CellPhoneDB output
  cpdb_pairs <- unique(immune_data$interacting_pair)
  cpdb_normalized <- sapply(cpdb_pairs, normalize_lr_pair)
  names(cpdb_normalized) <- cpdb_pairs

  # Match curated pairs to CellPhoneDB output (try both orderings)
  matched_pairs <- character(0)
  for (cp in curated_pairs) {
    cp_norm <- normalize_lr_pair(cp)
    # Find matching CPDB pair by normalized name
    matches <- names(cpdb_normalized)[cpdb_normalized == cp_norm]
    if (length(matches) > 0) {
      matched_pairs <- c(matched_pairs, matches[1])
    }
  }
  matched_pairs <- unique(matched_pairs)
  cat(sprintf("  Matched %d curated pairs via normalization\n", length(matched_pairs)))

  # Show unmatched curated pairs for debugging
  curated_normalized <- sapply(curated_pairs, normalize_lr_pair)
  unmatched_curated <- curated_pairs[!curated_normalized %in% cpdb_normalized]
  if (length(unmatched_curated) > 0 && length(unmatched_curated) <= 10) {
    cat(sprintf("  Unmatched curated pairs: %s\n", paste(unmatched_curated, collapse = ", ")))
  } else if (length(unmatched_curated) > 10) {
    cat(sprintf("  %d curated pairs not found in CellPhoneDB output\n", length(unmatched_curated)))
  }

  # Filter using matched pairs (not original curated_pairs)
  sig_lr_pairs <- immune_data %>%
    filter(interacting_pair %in% matched_pairs) %>%
    filter(!is.na(pvalue_fdr) & pvalue_fdr < 0.05) %>%
    group_by(interacting_pair) %>%
    summarise(
      n_sig = n(),
      n_young_sig = sum(age_group == "Young"),
      n_elderly_sig = sum(age_group == "Elderly"),
      min_pval = min(pvalue_fdr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_sig >= 1) %>%
    arrange(min_pval)

  cat(sprintf("  Matched curated L-R pairs: %d\n", nrow(sig_lr_pairs)))

  # Also show which curated pairs were NOT found (for debugging)
  missing_pairs <- setdiff(curated_pairs, sig_lr_pairs$interacting_pair)
  if (length(missing_pairs) > 0) {
    cat(sprintf("  Note: %d curated pairs not significant in Macrophage-immune data\n", length(missing_pairs)))
  }

} else {
  # Default: top50 mode - statistical filtering with FDR correction
  sig_lr_pairs <- immune_data %>%
    filter(!is.na(pvalue_fdr) & pvalue_fdr < 0.05) %>%
    group_by(interacting_pair) %>%
    summarise(
      n_sig = n(),
      n_young_sig = sum(age_group == "Young"),
      n_elderly_sig = sum(age_group == "Elderly"),
      min_pval = min(pvalue_fdr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # For macrophage-focused analysis, include pairs significant in either group
    filter(n_sig >= 1) %>%
    arrange(min_pval) %>%
    head(40)  # Top 40 L-R pairs (match manuscript ~40 rows)
}

# NOTE: Manuscript used a different CellPhoneDB database version with different L-R pair names
# Our regenerated figure uses current CellPhoneDB v5 naming conventions
# The biological patterns (Macrophage-immune communication by age) should be comparable

cat(sprintf("  Selected L-R pairs: %d\n", nrow(sig_lr_pairs)))

# Use ALL Macrophage|target cell pairs (manuscript shows all 12: 6 targets × 2 age groups)
sig_cell_pairs <- immune_data %>%
  filter(!is.na(pvalue),
         interacting_pair %in% sig_lr_pairs$interacting_pair) %>%
  distinct(cell_pair) %>%
  pull(cell_pair)

cat(sprintf("  Cell pairs: %d\n", length(sig_cell_pairs)))

# Create final plot data
plot_data <- immune_data %>%
  filter(interacting_pair %in% sig_lr_pairs$interacting_pair,
         cell_pair %in% sig_cell_pairs) %>%
  mutate(
    neg_log_pval = -log10(pvalue_fdr + 1e-10),
    neg_log_pval = pmin(neg_log_pval, 3),  # Cap at 3 to match manuscript scale
    # Use log2 mean expression to match manuscript
    log2_mean = log2(mean_expr + 1),
    # Set non-significant (FDR >= 0.05) to NA for cleaner plot
    log2_mean_plot = ifelse(pvalue_fdr < 0.05, log2_mean, NA)
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

# FIX: Order cell pairs by age using explicit factor levels (not alphabetical)
# "Younger" should appear BEFORE "Older" for each cell pair
plot_data <- plot_data %>%
  mutate(age_label = factor(age_label, levels = c("Younger", "Older"))) %>%
  arrange(cell_pair, age_label)

# Create ordered factor for cell_pair_age (preserves Younger-before-Older ordering)
cell_pair_age_order <- unique(plot_data$cell_pair_age)
plot_data$cell_pair_age <- factor(plot_data$cell_pair_age, levels = cell_pair_age_order)

# Create CellPhoneDB-style dot plot matching manuscript Figure 7D EXACTLY
# Key features: dense layout, top x-axis labels, blue-yellow-red gradient, SOLID dots
p <- ggplot(plot_data, aes(x = cell_pair_age, y = interacting_pair)) +
  geom_point(aes(size = neg_log_pval, color = log2_mean_plot), na.rm = TRUE) +  # SOLID dots, no outline
  scale_size_continuous(
    name = expression(-log[10](p-value)),
    range = c(2, 8),  # Larger dots like manuscript
    breaks = c(0, 1, 2, 3),
    limits = c(0, 3),
    guide = guide_legend(order = 1)
  ) +
  scale_color_gradientn(
    name = expression(Log[2]~mean~(Molecule~1~","~Molecule~2)),
    colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),  # RdBu diverging
    na.value = "grey90",
    limits = c(-10, 5),  # Manuscript range: -10 to +5
    guide = guide_colorbar(order = 2, barheight = unit(4, "cm"))
  ) +
  scale_x_discrete(position = "top") +  # X-axis labels on TOP (manuscript style)
  theme_minimal(base_size = 10) +
  theme(
    # X-axis: rotated labels on top - LARGER font
    axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 9, color = "black"),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    # Y-axis: L-R pair names - LARGER font
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_blank(),
    # Grid: very subtle
    panel.grid.major = element_line(color = "grey90", linewidth = 0.15),
    panel.grid.minor = element_blank(),
    # Panel border (manuscript has box around plot area)
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # Background
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # Legend: right side
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm"),
    # Tight margins for dense look
    plot.margin = margin(5, 5, 5, 5, "pt")
  )

# Save with mode-specific filename - PORTRAIT orientation like manuscript
mode_suffix <- ifelse(filter_mode == "curated", "_curated", "")
output_file <- file.path(figures_dir, sprintf("fig7d_cellphonedb_dotplot%s.png", mode_suffix))
ggsave(output_file, p, width = 8, height = 14, dpi = 300)  # Portrait: taller than wide
cat(sprintf("  Saved: %s\n", output_file))

# Also save as PDF
pdf_file <- file.path(figures_dir, sprintf("fig7d_cellphonedb_dotplot%s.pdf", mode_suffix))
ggsave(pdf_file, p, width = 8, height = 14)
cat(sprintf("  Saved: %s\n", pdf_file))

# Also save as SVG
svg_file <- file.path(figures_dir, sprintf("fig7d_cellphonedb_dotplot%s.svg", mode_suffix))
ggsave(svg_file, p, width = 8, height = 14, device = "svg")
cat(sprintf("  Saved: %s\n", svg_file))

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

output_csv <- file.path(script_dir, "outputs", sprintf("cellphonedb_summary%s.csv", mode_suffix))
fwrite(summary_df, output_csv)
cat(sprintf("  Saved: %s\n", output_csv))

cat(sprintf("\n=== CellPhoneDB comparative dot plot complete (mode: %s) ===\n", filter_mode))
