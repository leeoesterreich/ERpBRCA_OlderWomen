#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/08_run_wcsea.R
# WCSEA pathway analysis for single-cell RNA-seq
# Using indepthPathway package (https://github.com/wangxlab/indepthPathway)
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/wcsea_results.rds
#   - analysis/04_human_scrnaseq/outputs/wcsea_enrichment.csv

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

# Check if indepthPathway is installed
indepthPathway_available <- FALSE
if (requireNamespace("indepthPathway", quietly = TRUE)) {
  library(indepthPathway)
  indepthPathway_available <- TRUE
} else {
  cat("indepthPathway not available - will save markers for manual analysis\n")
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
output_dir <- file.path(script_dir, "outputs")

cat("=== WCSEA Pathway Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))

# Get normalized expression matrix
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
expr_matrix <- GetAssayData(seurat_obj, slot = "data")

cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Identify marker genes for comparison
# -----------------------------------------------------------------------------
cat("\nStep 2: Finding markers (Young vs Elderly)...\n")

# Subset to Young and Elderly only
seurat_subset <- subset(seurat_obj, AgeGroup %in% c("Young", "Elderly"))
Idents(seurat_subset) <- seurat_subset$AgeGroup

# Find differentially expressed genes
markers <- FindMarkers(
  seurat_subset,
  ident.1 = "Elderly",
  ident.2 = "Young",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# Apply FDR correction (already done by FindMarkers via p_val_adj)
markers_sig <- markers %>%
  filter(p_val_adj < 0.05) %>%
  tibble::rownames_to_column("gene")

cat("  Significant markers:", nrow(markers_sig), "\n")
cat("    Upregulated in Elderly:", sum(markers_sig$avg_log2FC > 0), "\n")
cat("    Downregulated in Elderly:", sum(markers_sig$avg_log2FC < 0), "\n")

# -----------------------------------------------------------------------------
# Step 3: Run WCSEA (if available)
# -----------------------------------------------------------------------------
cat("\nStep 3: Running WCSEA...\n")

# Prepare gene list with fold changes
gene_fc <- setNames(markers_sig$avg_log2FC, markers_sig$gene)

if (indepthPathway_available) {
  # Run WCSEA (Weighted Concept Signature Enrichment Analysis)
  # This is designed to handle noisy single-cell data
  tryCatch({
    wcsea_result <- WCSEA(
      gene_list = gene_fc,
      species = "human",
      pathway_db = "GO_BP",  # Gene Ontology Biological Process
      min_genes = 10,
      max_genes = 500
    )

    # Extract results
    wcsea_df <- wcsea_result$enrichment_result

    # FDR correction
    wcsea_df$padj <- p.adjust(wcsea_df$pvalue, method = "BH")
    wcsea_sig <- wcsea_df %>%
      filter(padj < 0.05) %>%
      arrange(padj)

    cat("  Enriched pathways (FDR < 0.05):", nrow(wcsea_sig), "\n")

    # Save results
    saveRDS(wcsea_result, file.path(output_dir, "wcsea_results.rds"))
    write.csv(wcsea_sig, file.path(output_dir, "wcsea_enrichment.csv"), row.names = FALSE)

  }, error = function(e) {
    cat("  Error running WCSEA:", conditionMessage(e), "\n")
    cat("  Saving marker genes for manual WCSEA analysis\n")
    write.csv(markers_sig, file.path(output_dir, "wcsea_input_markers.csv"), row.names = FALSE)
  })
} else {
  cat("  indepthPathway not available - saving marker genes instead\n")
  write.csv(markers_sig, file.path(output_dir, "wcsea_input_markers.csv"), row.names = FALSE)
}

# -----------------------------------------------------------------------------
# Step 4: Also check estrogen-related pathways specifically
# -----------------------------------------------------------------------------
cat("\nStep 4: Checking estrogen pathway genes...\n")

estrogen_genes <- c("ESR1", "ESR2", "PGR", "GREB1", "TFF1", "AREG", "XBP1")
estrogen_in_markers <- markers_sig %>%
  filter(gene %in% estrogen_genes)

if (nrow(estrogen_in_markers) > 0) {
  cat("  Estrogen-related genes in markers:\n")
  print(estrogen_in_markers[, c("gene", "avg_log2FC", "p_val_adj")])
} else {
  cat("  No canonical estrogen genes in significant markers\n")
}

cat("\n=== WCSEA complete ===\n")
