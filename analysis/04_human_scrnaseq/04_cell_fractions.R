#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/04_cell_fractions.R
# Cell type fraction analysis across age groups
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#   - analysis/04_human_scrnaseq/outputs/clinical_data.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/cell_fractions.csv
#   - analysis/04_human_scrnaseq/outputs/cell_fraction_stats.csv
#   - analysis/04_human_scrnaseq/outputs/fraction_boxplot.pdf
#   - analysis/04_human_scrnaseq/figures/fraction_boxplot.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
})

# Set Arial as default font for all plots
library(showtext)
# Use system Arial if available; showtext_auto() enables it for all devices
tryCatch(
  font_add("Arial", "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"),
  error = function(e) message("Arial not found; using default sans font")
)
showtext_auto()
theme_set(theme_bw(base_size = 14, base_family = "sans"))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "outputs")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Cell Fraction Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))
clinical_data <- readRDS(file.path(output_dir, "clinical_data.rds"))

cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Calculate cell type fractions per patient
# -----------------------------------------------------------------------------
cat("\nStep 2: Calculating cell fractions...\n")

# Count cells per patient per cell type
cell_counts <- table(seurat_obj$CellTypeAnnot, seurat_obj$orig.ident) %>%
  as.data.frame.matrix()

# Calculate proportions
cell_fractions <- sweep(cell_counts, 2, colSums(cell_counts), "/") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("CaseID")

# Merge with clinical data
fraction_data <- inner_join(
  clinical_data[, c("CaseID", "Age", "AgeGroup")],
  cell_fractions,
  by = "CaseID"
)

write.csv(fraction_data, file.path(output_dir, "cell_fractions.csv"), row.names = FALSE)
cat("  Saved: cell_fractions.csv\n")

# -----------------------------------------------------------------------------
# Step 3: Statistical tests with FDR correction
# -----------------------------------------------------------------------------
cat("\nStep 3: Running statistical tests...\n")

cell_types <- colnames(cell_fractions)[-1]  # Exclude CaseID

# Compare Young vs Elderly for each cell type
# FIX: Apply FDR correction across all comparisons
stat_results <- lapply(cell_types, function(ct) {
  young_vals <- fraction_data %>% filter(AgeGroup == "Young") %>% pull(!!sym(ct))
  elderly_vals <- fraction_data %>% filter(AgeGroup == "Elderly") %>% pull(!!sym(ct))

  # Wilcoxon test (non-parametric, appropriate for small n)
  if (length(young_vals) >= 2 && length(elderly_vals) >= 2) {
    test_result <- wilcox.test(young_vals, elderly_vals, exact = FALSE)
    data.frame(
      CellType = ct,
      Young_mean = mean(young_vals),
      Elderly_mean = mean(elderly_vals),
      pvalue = test_result$p.value
    )
  } else {
    data.frame(
      CellType = ct,
      Young_mean = mean(young_vals),
      Elderly_mean = mean(elderly_vals),
      pvalue = NA
    )
  }
}) %>% bind_rows()

# BIOSTATISTICAL FIX: Apply BH-FDR correction
stat_results$padj <- p.adjust(stat_results$pvalue, method = "BH")
stat_results$significant <- stat_results$padj < 0.05

write.csv(stat_results, file.path(output_dir, "cell_fraction_stats.csv"), row.names = FALSE)

cat("  Cell types tested:", nrow(stat_results), "\n")
cat("  Significant (FDR < 0.05):", sum(stat_results$significant, na.rm = TRUE), "\n")

# -----------------------------------------------------------------------------
# Step 4: Visualization
# -----------------------------------------------------------------------------
cat("\nStep 4: Generating plots...\n")

# Reshape for plotting
plot_data <- fraction_data %>%
  pivot_longer(
    cols = all_of(cell_types),
    names_to = "CellType",
    values_to = "Fraction"
  ) %>%
  filter(AgeGroup != "MidAge")  # Compare Young vs Elderly

p <- ggplot(plot_data, aes(x = CellType, y = Fraction, fill = AgeGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = AgeGroup), position = position_jitterdodge(jitter.width = 0.1)) +
  scale_fill_manual(values = c("Young" = "#56B4E9", "Elderly" = "#D55E00")) +
  scale_color_manual(values = c("Young" = "#56B4E9", "Elderly" = "#D55E00")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(x = "", y = "Cell Type Fraction", title = "Cell Type Proportions by Age Group")

pdf(file.path(output_dir, "fraction_boxplot.pdf"), width = 12, height = 6)
print(p)
dev.off()

# Save PNG for validation pipeline
ggsave(file.path(figures_dir, "fraction_boxplot.png"), p,
       width = 12, height = 6, dpi = 300, bg = "white")
ggsave(file.path(figures_dir, "fraction_boxplot.svg"), p,
       width = 12, height = 6)  # SVG for vector assembly

cat("\n=== Cell fraction analysis complete ===\n")
