#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/08_differential_expression.R
# Differential expression analysis: Young vs Aged per cell type
# BIOSTATISTICAL FIX: This analysis was missing from original
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
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
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

cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  Age groups:\n")
print(table(seurat_obj$AgeGroup))

# Ensure AgeGroup is set correctly
if (is.null(seurat_obj$AgeGroup) || all(is.na(seurat_obj$AgeGroup))) {
  stop("ERROR: AgeGroup metadata not found or all NA")
}

# Set RNA assay as default for DE
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

# -----------------------------------------------------------------------------
# Step 2: Run DE for Each Cell Type
# -----------------------------------------------------------------------------
cat("\nStep 2: Running differential expression per cell type...\n")

cell_types <- unique(seurat_obj$CellTypeByMarker_RatsnRNAseq)
cell_types <- cell_types[!is.na(cell_types)]
all_de_results <- data.frame()

cat("  Cell types to analyze:", length(cell_types), "\n")
cat("    ", paste(cell_types, collapse = ", "), "\n\n")

for (ct in cell_types) {
  cat("  Processing:", ct, "\n")

  # Subset to cell type
  seurat_ct <- subset(seurat_obj, CellTypeByMarker_RatsnRNAseq == ct)

  # Check sample sizes
  n_aged <- sum(seurat_ct$AgeGroup == "Aged")
  n_young <- sum(seurat_ct$AgeGroup == "Young")
  cat("    Aged:", n_aged, "| Young:", n_young, "\n")

  # Skip if too few cells
  if (n_aged < 10 || n_young < 10) {
    cat("    Skipping - too few cells (min 10 per group required)\n")
    next
  }

  # Set identity to AgeGroup
  Idents(seurat_ct) <- "AgeGroup"

  # Run FindMarkers (Wilcoxon test with BH correction)
  tryCatch({
    de_results <- FindMarkers(
      seurat_ct,
      ident.1 = "Aged",
      ident.2 = "Young",
      test.use = "wilcox",
      min.pct = 0.1,
      logfc.threshold = 0.25,
      verbose = FALSE
    )

    if (nrow(de_results) > 0) {
      de_results$gene <- rownames(de_results)
      de_results$celltype <- ct
      # BIOSTATISTICAL FIX: Apply BH FDR correction
      de_results$FDR <- p.adjust(de_results$p_val, method = "BH")
      all_de_results <- rbind(all_de_results, de_results)
      cat("    DE genes:", nrow(de_results), "\n")
    } else {
      cat("    No DE genes found\n")
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
  arrange(FDR, p_val) %>%
  mutate(
    Sig_nominal = p_val < 0.05,
    Sig_FDR = FDR < 0.05,
    Direction = ifelse(avg_log2FC > 0, "Up_in_Aged", "Down_in_Aged")
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

pdf(file.path(fig_dir, "DE_volcano_plots.pdf"), width = 12, height = 10)

for (ct in unique(all_de_results$celltype)) {
  ct_results <- all_de_results %>% filter(celltype == ct)

  # Label top genes
  top_genes <- ct_results %>%
    filter(Sig_FDR) %>%
    slice_min(FDR, n = 10)

  p <- ggplot(ct_results, aes(x = avg_log2FC, y = -log10(FDR))) +
    geom_point(aes(color = Sig_FDR), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray") +
    geom_text_repel(data = top_genes, aes(label = gene), max.overlaps = 20, size = 3) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    theme_bw() +
    labs(
      title = paste("Differential Expression:", ct),
      subtitle = paste("Aged vs Young | FDR<0.05:", sum(ct_results$Sig_FDR)),
      x = "log2 Fold Change (Aged/Young)",
      y = "-log10(FDR)"
    ) +
    theme(legend.position = "none")

  print(p)
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

cat("  Volcano plots saved to:", file.path(fig_dir, "DE_volcano_plots.pdf"), "\n")

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
