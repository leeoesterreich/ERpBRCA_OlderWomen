#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/08_differential_expression.R
# Differential expression analysis: Young vs Aged per cell type
# Uses pseudobulk DESeq2 to avoid pseudoreplication (Squair et al. 2021)
#
# Inputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - results/corrected/rat_snrnaseq/DE_results_by_celltype.csv
#   - figures/by_analysis/rat_snrnaseq/DE_volcano_plots.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(future)
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

# Increase memory limit for parallelization (required for large cell populations)
options(future.globals.maxSize = 4 * 1024^3)  # 4 GB

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

cat("=== Differential Expression Analysis: Young vs Aged ===\n")
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
cat("Step 1: Loading annotated Seurat object...\n")
seurat_obj <- readRDS(input_rds)

# Validate Seurat object
if (!inherits(seurat_obj, "Seurat")) {
  stop("ERROR: Loaded object is not a valid Seurat object")
}
if (ncol(seurat_obj) == 0) {
  stop("ERROR: Seurat object contains no cells")
}

cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  Age groups:\n")
print(table(seurat_obj$AgeGroup))

# Ensure AgeGroup is set correctly
if (is.null(seurat_obj$AgeGroup) || all(is.na(seurat_obj$AgeGroup))) {
  stop("ERROR: AgeGroup metadata not found or all NA")
}

# Set RNA assay as default for DE
DefaultAssay(seurat_obj) <- "RNA"

# Seurat v5 requires joining layers before differential expression
if ("JoinLayers" %in% ls("package:SeuratObject")) {
  seurat_obj <- JoinLayers(seurat_obj)
  cat("  Joined Seurat v5 layers\n")
}

seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

# -----------------------------------------------------------------------------
# Step 2: Pseudobulk DE for Each Cell Type (DESeq2)
# -----------------------------------------------------------------------------
cat("\nStep 2: Running pseudobulk differential expression per cell type...\n")
cat("  Method: DESeq2 on aggregated counts per sample (avoids pseudoreplication)\n\n")

cell_types <- unique(seurat_obj$CellTypeByMarker_RatsnRNAseq)
cell_types <- cell_types[!is.na(cell_types)]
all_de_results <- data.frame()

MIN_CELLS_PER_SAMPLE <- 10
MIN_SAMPLES_PER_GROUP <- 2

cat("  Cell types to analyze:", length(cell_types), "\n")
cat("    ", paste(cell_types, collapse = ", "), "\n\n")

for (ct in cell_types) {
  cat("  Processing:", ct, "\n")

  # Subset to cell type
  seurat_ct <- subset(seurat_obj, CellTypeByMarker_RatsnRNAseq == ct)

  # Check cells per sample and filter
  cells_per_sample <- table(seurat_ct$orig.ident)
  valid_samples <- names(cells_per_sample[cells_per_sample >= MIN_CELLS_PER_SAMPLE])

  if (length(valid_samples) == 0) {
    cat("    Skipping - no samples with >=", MIN_CELLS_PER_SAMPLE, "cells\n")
    next
  }

  seurat_ct <- subset(seurat_ct, orig.ident %in% valid_samples)

  # Check biological replicates per group
  sample_groups <- unique(data.frame(
    sample = seurat_ct$orig.ident,
    group = seurat_ct$AgeGroup,
    stringsAsFactors = FALSE
  ))
  n_aged_samples <- sum(sample_groups$group == "Aged")
  n_young_samples <- sum(sample_groups$group == "Young")
  cat("    Aged samples:", n_aged_samples, "| Young samples:", n_young_samples, "\n")

  if (n_aged_samples < MIN_SAMPLES_PER_GROUP || n_young_samples < MIN_SAMPLES_PER_GROUP) {
    cat("    Skipping - need >=", MIN_SAMPLES_PER_GROUP, "samples per group\n")
    next
  }

  # Aggregate counts per sample (pseudobulk)
  tryCatch({
    pseudo_counts <- AggregateExpression(
      seurat_ct,
      group.by = "orig.ident",
      assays = "RNA",
      slot = "counts",
      return.seurat = FALSE
    )$RNA

    # Build sample metadata for DESeq2
    sample_meta <- sample_groups[!duplicated(sample_groups$sample), ]
    rownames(sample_meta) <- sample_meta$sample
    sample_meta <- sample_meta[colnames(pseudo_counts), , drop = FALSE]
    sample_meta$group <- factor(sample_meta$group, levels = c("Young", "Aged"))

    # Filter low-count genes: require >= 10 counts in >= 2 samples
    keep <- rowSums(pseudo_counts >= 10) >= 2
    pseudo_counts <- pseudo_counts[keep, ]

    if (nrow(pseudo_counts) < 10) {
      cat("    Skipping - too few genes pass filter (", nrow(pseudo_counts), ")\n")
      next
    }

    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(
      countData = pseudo_counts,
      colData = sample_meta,
      design = ~ group
    )
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("group", "Aged", "Young"), alpha = 0.05)

    de_results <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      filter(!is.na(padj)) %>%
      mutate(
        celltype = ct,
        FDR = padj,
        avg_log2FC = log2FoldChange
      )

    if (nrow(de_results) > 0) {
      all_de_results <- rbind(all_de_results, de_results)
      n_sig <- sum(de_results$FDR < 0.05, na.rm = TRUE)
      cat("    Tested genes:", nrow(de_results), "| Significant (FDR<0.05):", n_sig, "\n")
    } else {
      cat("    No genes passed filters\n")
    }
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
  })
}

# Check if we have any results
if (nrow(all_de_results) == 0) {
  warning("No DE results found. Check data quality and sample sizes.")
  cat("\nWARNING: No differential expression results generated.\n")
  cat("This may indicate:\n")
  cat("  - Too few cells per cell type\n")
  cat("  - Insufficient biological variation\n")
  cat("  - Data quality issues\n")
  quit(status = 0)
}

# -----------------------------------------------------------------------------
# Step 3: Summarize Results
# -----------------------------------------------------------------------------
cat("\nStep 3: Summarizing results...\n")

all_de_results <- all_de_results %>%
  arrange(FDR, pvalue) %>%
  mutate(
    Sig_nominal = pvalue < 0.05,
    Sig_FDR = FDR < 0.05,
    Direction = ifelse(avg_log2FC > 0, "Up_in_Aged", "Down_in_Aged"),
    # FC-filtered significance: require |log2FC|>0.5 in addition to FDR<0.05
    # Prevents cell-count artifacts (e.g. Luminal n=3 pseudobulk) from flooding plots
    Sig_FC = case_when(
      FDR < 0.05 & avg_log2FC > 0.5  ~ "Up",
      FDR < 0.05 & avg_log2FC < -0.5 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Summary by cell type
de_summary <- all_de_results %>%
  group_by(celltype) %>%
  summarise(
    total_DE_genes = n(),
    sig_nominal = sum(Sig_nominal),
    sig_FDR = sum(Sig_FDR),
    up_in_aged = sum(Sig_FDR & Direction == "Up_in_Aged"),
    down_in_aged = sum(Sig_FDR & Direction == "Down_in_Aged"),
    .groups = "drop"
  )

cat("\nDE Summary by Cell Type:\n")
print(as.data.frame(de_summary))

# -----------------------------------------------------------------------------
# Step 4: Generate Volcano Plots
# -----------------------------------------------------------------------------
cat("\nStep 4: Generating volcano plots...\n")

# Store individual plots for combined figure
volcano_plots <- list()

pdf(file.path(fig_dir, "DE_volcano_plots.pdf"), width = 12, height = 10)

for (ct in unique(all_de_results$celltype)) {
  ct_results <- all_de_results %>% filter(celltype == ct)

  # Label top genes (require FC filter so labels come from biologically meaningful hits)
  top_genes <- ct_results %>%
    filter(Sig_FC != "NS") %>%
    slice_min(FDR, n = 5)

  p <- ggplot(ct_results, aes(x = avg_log2FC, y = -log10(FDR))) +
    geom_point(aes(color = Sig_FC), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
    geom_text_repel(data = top_genes, aes(label = gene), max.overlaps = 10, size = 3.2) +
    scale_color_manual(values = c("NS" = "#BDBDBD", "Up" = "#D55E00", "Down" = "#0072B2")) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Differential Expression:", ct),
      subtitle = sprintf("%d up (|FC|>1.5), %d down at FDR<0.05",
        sum(ct_results$Sig_FC == "Up",   na.rm = TRUE),
        sum(ct_results$Sig_FC == "Down", na.rm = TRUE)),
      x = "log2 Fold Change (Aged/Young)",
      y = "-log10(FDR)"
    ) +
    theme(legend.position = "none")

  # Cap y-axis (extreme -log10(FDR) from high cell counts makes plot unreadable)
  y_cap <- min(50, max(-log10(ct_results$FDR + 1e-300), na.rm = TRUE))
  p <- p + ylim(0, y_cap)

  print(p)
  volcano_plots[[ct]] <- p

  # Save individual volcano plots as PNG and SVG
  ct_clean <- gsub("[^A-Za-z0-9]", "_", ct)
  ggsave(file.path(fig_dir, paste0("DE_volcano_", ct_clean, ".png")), p, width = 10, height = 8, dpi = 300)
  ggsave(file.path(fig_dir, paste0("DE_volcano_", ct_clean, ".svg")), p, width = 10, height = 8)
}

# Combined summary plot
p_summary <- ggplot(de_summary, aes(x = reorder(celltype, sig_FDR), y = sig_FDR)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sig_FDR), hjust = -0.2) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Significant DE Genes by Cell Type",
    subtitle = "FDR < 0.05",
    x = "Cell Type",
    y = "Number of DE Genes"
  )
print(p_summary)

dev.off()

# Save summary plot as PNG and SVG
ggsave(file.path(fig_dir, "DE_summary_barplot.png"), p_summary, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "DE_summary_barplot.svg"), p_summary, width = 8, height = 6)

# Create combined 2x2 volcano plot grid
if (length(volcano_plots) == 4) {
  library(patchwork)
  p_combined <- (volcano_plots[[1]] + volcano_plots[[2]]) / (volcano_plots[[3]] + volcano_plots[[4]]) +
    plot_annotation(title = "Differential Expression: Young vs Aged by Cell Type")
  ggsave(file.path(fig_dir, "DE_volcano_combined.png"), p_combined, width = 16, height = 14, dpi = 300)
  ggsave(file.path(fig_dir, "DE_volcano_combined.svg"), p_combined, width = 16, height = 14)
}

cat("  Volcano plots saved to:", file.path(fig_dir, "DE_volcano_plots.pdf"), "\n")
cat("  Individual PNG/SVG plots also saved\n")

# -----------------------------------------------------------------------------
# Step 5: Save Results
# -----------------------------------------------------------------------------
cat("\nStep 5: Saving results...\n")

write.csv(all_de_results, file.path(results_dir, "DE_results_by_celltype.csv"), row.names = FALSE)
write.csv(de_summary, file.path(results_dir, "DE_summary_by_celltype.csv"), row.names = FALSE)

# Save top DE genes per cell type
top_de_genes <- all_de_results %>%
  filter(Sig_FDR) %>%
  group_by(celltype) %>%
  slice_min(FDR, n = 20) %>%
  ungroup()

write.csv(top_de_genes, file.path(results_dir, "DE_top20_genes_by_celltype.csv"), row.names = FALSE)

cat("  Results saved:\n")
cat("    - DE_results_by_celltype.csv (all results)\n")
cat("    - DE_summary_by_celltype.csv (summary statistics)\n")
cat("    - DE_top20_genes_by_celltype.csv (top 20 per cell type)\n")

cat("\n=== Differential expression analysis complete ===\n")
cat("Results:", file.path(results_dir, "DE_results_by_celltype.csv"), "\n")
