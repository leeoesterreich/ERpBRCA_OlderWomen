#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/05_gene_expression_violin.R
# Gene expression violin plots by age group and cell type
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/violin_*.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
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
output_dir <- file.path(script_dir, "outputs")

cat("=== Gene Expression Violin Plots ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load and prepare data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))

# Normalize if not already done
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

# Order age groups
seurat_obj$AgeGroup <- factor(seurat_obj$AgeGroup, levels = c("Young", "MidAge", "Elderly"))

cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Define gene sets
# -----------------------------------------------------------------------------
cat("\nStep 2: Defining gene sets...\n")

immune_checkpoint <- c("CTLA4", "PDCD1", "CD274", "LAG3", "HAVCR2", "PDCD1LG2")
macrophage_markers <- c("CCL2", "CCL3", "CCL4", "TNF", "TGFB1", "CD163", "MRC1")
cytokines <- c("IL1A", "IL1B", "IL6", "CXCL8", "CXCL10", "CCL5")

# Filter to genes present in data
filter_genes <- function(genes) {
  present <- genes[genes %in% rownames(seurat_obj)]
  if (length(present) < length(genes)) {
    missing <- setdiff(genes, present)
    cat("  Missing:", paste(missing, collapse = ", "), "\n")
  }
  present
}

immune_checkpoint <- filter_genes(immune_checkpoint)
macrophage_markers <- filter_genes(macrophage_markers)
cytokines <- filter_genes(cytokines)

# -----------------------------------------------------------------------------
# Step 3: Generate violin plots by cell type
# -----------------------------------------------------------------------------
cat("\nStep 3: Generating plots...\n")

cell_types <- unique(seurat_obj$CellTypeAnnot)

for (ct in cell_types) {
  cat("  Processing:", ct, "\n")

  # Subset to cell type
  seurat_subset <- subset(seurat_obj, CellTypeAnnot == ct)

  if (ncol(seurat_subset) < 10) {
    cat("    Skipping (too few cells)\n")
    next
  }

  # Immune checkpoint genes
  if (length(immune_checkpoint) > 0) {
    p <- VlnPlot(
      seurat_subset,
      features = immune_checkpoint,
      group.by = "AgeGroup",
      pt.size = 0
    ) + plot_annotation(title = paste(ct, "- Immune Checkpoint"))

    pdf(file.path(output_dir, paste0("violin_checkpoint_", ct, ".pdf")), width = 12, height = 8)
    print(p)
    dev.off()
  }

  # Macrophage markers (only for macrophage/monocyte)
  if (grepl("Macro|Mono", ct, ignore.case = TRUE) && length(macrophage_markers) > 0) {
    p <- VlnPlot(
      seurat_subset,
      features = macrophage_markers,
      group.by = "AgeGroup",
      pt.size = 0
    ) + plot_annotation(title = paste(ct, "- Macrophage Markers"))

    pdf(file.path(output_dir, paste0("violin_macrophage_", ct, ".pdf")), width = 12, height = 8)
    print(p)
    dev.off()
  }
}

cat("\n=== Violin plots complete ===\n")
