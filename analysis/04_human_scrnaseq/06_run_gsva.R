#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/06_run_gsva.R
# GSVA on scRNA-seq data — supports both pseudo-bulk and single-cell modes
#
# Manuscript methods (line 226-227):
#   "GSVA (version 1.48.3) with default parameters"
#   "gene set collections associated with the estrogen pathway ... sourced from
#    hallmark, reactome, wikipathways, and gene ontology biological pathways"
#
# Usage:
#   Rscript 06_run_gsva.R                    # runs BOTH modes
#   Rscript 06_run_gsva.R --mode=pseudobulk  # pseudo-bulk only
#   Rscript 06_run_gsva.R --mode=singlecell  # single-cell only (original approach)
#
# The original code (Sanghoon's GitHub, Step 6 second instance) does:
#   1. rowSums(counts) > 10  (gene filter)
#   2. Remove non-protein-coding genes (containing "." in gene name)
#   3. colSums > 1000 (cell filter)
#   → 18,063 genes x 28,732 cells
#   4. Run GSVA with default parameters
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/gsva_pseudobulk.rds
#   - analysis/04_human_scrnaseq/outputs/gsva_singlecell.rds
#   - analysis/04_human_scrnaseq/outputs/gsva_heatmap.pdf
#   - analysis/04_human_scrnaseq/outputs/gsva_comparison.csv
#   - analysis/04_human_scrnaseq/figures/gsva_heatmap*.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(GSVA)
  library(msigdbr)
  library(dplyr)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(pheatmap)
  library(tibble)
})

# Parse --mode argument
args <- commandArgs(trailingOnly = TRUE)
mode_arg <- grep("--mode=", args, value = TRUE)
if (length(mode_arg) > 0) {
  run_mode <- sub("--mode=", "", mode_arg)
} else {
  run_mode <- "both"
}
cat("=== scRNA-seq GSVA (mode:", run_mode, ") ===\n")

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
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))
cat("  Cells:", ncol(seurat_obj), "\n")

DefaultAssay(seurat_obj) <- "RNA"
counts_all <- GetAssayData(seurat_obj, slot = "counts")

# -----------------------------------------------------------------------------
# Step 2: Load gene sets
# Manuscript: "hallmark, reactome, wikipathways, and gene ontology biological pathways"
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading gene sets...\n")

all_sets <- msigdbr(species = "Homo sapiens")
all_list <- split(all_sets$gene_symbol, all_sets$gs_name)

find_pathway <- function(pattern, verbose = TRUE) {
  hit <- grep(pattern, names(all_list), ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) {
    if (verbose) cat("  WARNING: No match for pattern '", pattern, "'\n", sep = "")
    return(NULL)
  }
  if (verbose) cat("  Found: ", hit[1], " (", length(all_list[[hit[1]]]), " genes)\n", sep = "")
  all_list[[hit[1]]]
}

cat("  Looking up estrogen pathways...\n")
estrogen_pathways <- list(
  # WikiPathways
  "Estrogen metabolism WP697" = find_pathway("WP_ESTROGEN_METABOLISM_WP697$"),
  "Estrogen metabolism WP5276" = find_pathway("WP_ESTROGEN_METABOLISM_WP5276$"),
  "Estrogen signaling WP" = find_pathway("WP_ESTROGEN_SIGNALING"),
  # HALLMARK
  "Estrogen response late" = find_pathway("HALLMARK_ESTROGEN_RESPONSE_LATE"),
  "Estrogen response early" = find_pathway("HALLMARK_ESTROGEN_RESPONSE_EARLY"),
  # Reactome (all 6 estrogen sets)
  "Estrogen biosynthesis" = find_pathway("REACTOME_ESTROGEN_BIOSYNTHESIS$"),
  "Estrogen dependent gene expression" = find_pathway("REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION$"),
  "Estrogen dependent nuclear events" = find_pathway("REACTOME_ESTROGEN_DEPENDENT_NUCLEAR_EVENTS"),
  "Estrogen stimulated signaling PRKCZ" = find_pathway("REACTOME_ESTROGEN_STIMULATED_SIGNALING"),
  "Extra nuclear estrogen signaling" = find_pathway("REACTOME_EXTRA_NUCLEAR_ESTROGEN"),
  "RUNX1 regulates ER transcription" = find_pathway("REACTOME_RUNX1_REGULATES_ESTROGEN"),
  # GO Biological Process
  "Response to estrogen" = find_pathway("GOBP_RESPONSE_TO_ESTROGEN$"),
  "Estrogen receptor signaling pathway" = find_pathway("GOBP_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY"),
  "Positive regulation of ER signaling" = find_pathway("GOBP_POSITIVE_REGULATION_OF_INTRACELLULAR_ESTROGEN"),
  "Cellular response to estrogen stimulus" = find_pathway("GOBP_CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS")
)
estrogen_pathways <- estrogen_pathways[!sapply(estrogen_pathways, is.null)]

# Add LI_ESTROGENE from local GMT files
gmt_dir <- file.path(project_root, "data", "gmt")
parse_gmt <- function(path) {
  line <- readLines(path, n = 1)
  fields <- strsplit(line, "\t")[[1]]
  fields[-(1:2)]
}
if (file.exists(file.path(gmt_dir, "LI_ESTROGENE_EARLY_E2_RESPONSE_UP.v2025.1.Hs.gmt"))) {
  estrogen_pathways[["LI EstroGene early E2 response up"]] <- parse_gmt(file.path(gmt_dir, "LI_ESTROGENE_EARLY_E2_RESPONSE_UP.v2025.1.Hs.gmt"))
  estrogen_pathways[["LI EstroGene late E2 response up"]] <- parse_gmt(file.path(gmt_dir, "LI_ESTROGENE_LATE_E2_RESPONSE_UP.v2025.1.Hs.gmt"))
}
cat("  Gene sets loaded:", length(estrogen_pathways), "\n")

# =============================================================================
# Helper: aggregate GSVA by HSD17B7 status and make heatmap
# =============================================================================
make_hsd17b7_heatmap <- function(gsva_mat, expr_mat, suffix, seurat_obj) {
  # Get HSD17B7 expression per sample
  hsd17b7_expr <- expr_mat["HSD17B7", colnames(gsva_mat)]
  hsd17b7_median <- median(hsd17b7_expr, na.rm = TRUE)
  hsd17b7_group <- ifelse(hsd17b7_expr > hsd17b7_median, "HSD17B7+", "HSD17B7-")

  cat("  HSD17B7 median:", round(hsd17b7_median, 2), "\n")
  cat("  HSD17B7- samples:", sum(hsd17b7_group == "HSD17B7-"), "\n")
  cat("  HSD17B7+ samples:", sum(hsd17b7_group == "HSD17B7+"), "\n")

  # Aggregate GSVA scores
  gsva_agg <- sapply(c("HSD17B7-", "HSD17B7+"), function(g) {
    cols <- names(hsd17b7_group)[hsd17b7_group == g]
    rowMeans(gsva_mat[, cols, drop = FALSE], na.rm = TRUE)
  })

  gsva_agg <- pmax(pmin(gsva_agg, 4), -4)
  valid_rows <- intersect(names(estrogen_pathways), rownames(gsva_agg))
  gsva_agg <- gsva_agg[valid_rows, , drop = FALSE]

  data_range <- max(abs(gsva_agg), na.rm = TRUE)
  if (data_range < 0.01) data_range <- 0.5  # avoid degenerate color scale

  col_fun <- colorRamp2(
    c(-data_range, -data_range/2, 0, data_range/2, data_range),
    c("#313695", "#74add1", "white", "#fdae61", "#d73027")
  )

  ht <- Heatmap(
    gsva_agg,
    name = "Pathway\nactivity",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 12),
    column_names_side = "bottom",
    rect_gp = gpar(col = NA),
    width = unit(2, "cm"),
    height = unit(nrow(gsva_agg) * 0.5, "cm"),
    heatmap_legend_param = list(
      title = "Pathway\nactivity",
      at = c(-round(data_range, 1), 0, round(data_range, 1)),
      labels = c(as.character(-round(data_range, 1)), "0", as.character(round(data_range, 1))),
      legend_height = unit(3, "cm")
    )
  )

  # Save heatmap
  pdf(file.path(output_dir, paste0("gsva_heatmap_", suffix, ".pdf")), width = 10, height = 6)
  draw(ht, padding = unit(c(2, 2, 2, 2), "cm"))
  dev.off()

  png(file.path(figures_dir, paste0("gsva_heatmap_", suffix, ".png")), width = 10*300, height = 6*300, res = 300)
  draw(ht, padding = unit(c(2, 2, 2, 2), "cm"))
  dev.off()

  # Per-sample heatmap
  age_groups <- seurat_obj@meta.data %>%
    select(orig.ident, AgeGroup) %>%
    distinct()
  rownames(age_groups) <- age_groups$orig.ident

  sample_cols <- intersect(colnames(gsva_mat), age_groups$orig.ident)
  if (length(sample_cols) > 1) {
    ann_col <- data.frame(
      AgeGroup = age_groups[sample_cols, "AgeGroup"],
      HSD17B7 = hsd17b7_group[sample_cols],
      row.names = sample_cols
    )
    ann_colors <- list(
      AgeGroup = c(Young = "#4DAF4A", MidAge = "#377EB8", Elderly = "#E41A1C"),
      HSD17B7 = c("HSD17B7-" = "#74add1", "HSD17B7+" = "#d73027")
    )

    png(file.path(figures_dir, paste0("gsva_heatmap_per_sample_", suffix, ".png")),
        width = 12*300, height = 8*300, res = 300)
    pheatmap(
      gsva_mat[, sample_cols],
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      cluster_rows = FALSE,
      # GSVA scores are already bounded [-1, 1]; do not apply additional row-scaling
      scale = "none",
      show_colnames = TRUE,
      main = paste0("GSVA per-sample (", suffix, ")")
    )
    dev.off()
  }

  return(list(gsva_mat = gsva_mat, gsva_agg = gsva_agg, hsd17b7_group = hsd17b7_group))
}

# =============================================================================
# MODE 1: Pseudo-bulk GSVA (refactored approach)
# =============================================================================
pseudobulk_result <- NULL
if (run_mode %in% c("both", "pseudobulk")) {
  cat("\n=== PSEUDO-BULK GSVA ===\n")

  # Aggregate counts per patient
  patients <- unique(seurat_obj$orig.ident)
  pseudobulk <- sapply(patients, function(pt) {
    cells <- colnames(seurat_obj)[seurat_obj$orig.ident == pt]
    rowSums(counts_all[, cells, drop = FALSE])
  })

  # CPM + log2
  pseudobulk_cpm <- sweep(pseudobulk, 2, colSums(pseudobulk), "/") * 1e6
  pseudobulk_log <- log2(pseudobulk_cpm + 1)

  cat("  Pseudo-bulk matrix (raw):", nrow(pseudobulk_log), "genes x", ncol(pseudobulk_log), "samples\n")

  # Gene filtering — match single-cell path: remove low-count genes and non-coding
  gene_sums <- rowSums(pseudobulk)
  pseudobulk_log <- pseudobulk_log[gene_sums > 10, ]
  cat("  After gene filter (rowSums > 10):", nrow(pseudobulk_log), "genes\n")

  # Remove non-protein-coding genes (containing "." in name, e.g., RP11-34P13.7)
  non_coding <- grepl("\\.", rownames(pseudobulk_log))
  pseudobulk_log <- pseudobulk_log[!non_coding, ]
  cat("  After removing non-protein-coding:", nrow(pseudobulk_log), "genes\n")

  # Verify gene set overlap before GSVA
  all_gs_genes <- unique(unlist(estrogen_pathways))
  overlap <- sum(all_gs_genes %in% rownames(pseudobulk_log))
  cat("  Gene set overlap:", overlap, "of", length(all_gs_genes), "genes\n")

  # GSVA with Gaussian kcdf (appropriate for log-transformed continuous data)
  cat("  Running GSVA (kcdf = Gaussian)...\n")
  gsva_pb <- gsva(
    gsvaParam(
      as.matrix(pseudobulk_log),
      estrogen_pathways,
      kcdf = "Gaussian",
      maxDiff = TRUE
    )
  )
  cat("  Result:", nrow(gsva_pb), "pathways x", ncol(gsva_pb), "samples\n")

  pseudobulk_result <- make_hsd17b7_heatmap(gsva_pb, pseudobulk_log, "pseudobulk", seurat_obj)

  saveRDS(gsva_pb, file.path(output_dir, "gsva_pseudobulk.rds"))
  saveRDS(pseudobulk_log, file.path(output_dir, "pseudobulk_log2cpm.rds"))
  write.csv(pseudobulk_result$gsva_agg, file.path(output_dir, "gsva_aggregated_pseudobulk.csv"))
}

# =============================================================================
# MODE 2: Single-cell GSVA (original approach, matching Sanghoon's code)
# =============================================================================
singlecell_result <- NULL
if (run_mode %in% c("both", "singlecell")) {
  cat("\n=== SINGLE-CELL GSVA (original approach) ===\n")

  # Step A: Gene filtering — rowSums > 10 (matches original)
  counts_sc <- as.matrix(counts_all)  # dense for filtering
  gene_counts <- rowSums(counts_sc)
  counts_sc <- counts_sc[gene_counts > 10, ]
  cat("  After gene filter (rowSums > 10):", nrow(counts_sc), "genes\n")

  # Step B: Remove non-protein-coding genes (containing "." in name)
  # Original: dplyr::filter(!grepl("\\.", rownames(...)))
  non_coding <- grepl("\\.", rownames(counts_sc))
  counts_sc <- counts_sc[!non_coding, ]
  cat("  After removing non-protein-coding:", nrow(counts_sc), "genes\n")

  # Step C: Cell filtering — colSums > 1000 (matches original)
  cell_counts <- colSums(counts_sc)
  counts_sc <- counts_sc[, cell_counts > 1000]
  cat("  After cell filter (colSums > 1000):", ncol(counts_sc), "cells\n")
  cat("  (Manuscript/original expects: ~18,063 genes x ~28,732 cells)\n")

  # Step D: Run GSVA with default parameters (manuscript: "GSVA with default parameters")
  # For raw counts, use Poisson kcdf (GSVA default for integer count data)
  cat("  Running GSVA (kcdf = Poisson, default params)...\n")
  gsva_sc <- gsva(
    gsvaParam(
      counts_sc,
      estrogen_pathways,
      kcdf = "Poisson",
      maxDiff = TRUE
    )
  )

  # gsva_sc is now pathways x cells — aggregate per patient for comparison
  cat("  Single-cell result:", nrow(gsva_sc), "pathways x", ncol(gsva_sc), "cells\n")

  # Aggregate per patient (mean GSVA score across cells)
  kept_cells <- colnames(gsva_sc)
  # Map cells back to patients
  cell_patient <- seurat_obj$orig.ident[kept_cells]
  patients <- unique(cell_patient)

  gsva_sc_by_patient <- sapply(patients, function(pt) {
    pt_cells <- kept_cells[cell_patient == pt]
    rowMeans(gsva_sc[, pt_cells, drop = FALSE], na.rm = TRUE)
  })

  cat("  Aggregated to:", nrow(gsva_sc_by_patient), "pathways x",
      ncol(gsva_sc_by_patient), "patients\n")

  # Need expression per patient for HSD17B7 — use aggregated counts
  pseudobulk_sc <- sapply(patients, function(pt) {
    pt_cells <- kept_cells[cell_patient == pt]
    if (length(pt_cells) == 1) {
      counts_sc[, pt_cells]
    } else {
      rowSums(counts_sc[, pt_cells, drop = FALSE])
    }
  })
  pseudobulk_sc_cpm <- sweep(pseudobulk_sc, 2, colSums(pseudobulk_sc), "/") * 1e6
  pseudobulk_sc_log <- log2(pseudobulk_sc_cpm + 1)

  singlecell_result <- make_hsd17b7_heatmap(gsva_sc_by_patient, pseudobulk_sc_log,
                                              "singlecell", seurat_obj)

  saveRDS(gsva_sc_by_patient, file.path(output_dir, "gsva_singlecell.rds"))
  write.csv(singlecell_result$gsva_agg, file.path(output_dir, "gsva_aggregated_singlecell.csv"))
}

# =============================================================================
# Comparison: if both modes ran, output side-by-side comparison
# =============================================================================
if (!is.null(pseudobulk_result) && !is.null(singlecell_result)) {
  cat("\n=== COMPARING PSEUDO-BULK vs SINGLE-CELL GSVA ===\n")

  pb_agg <- pseudobulk_result$gsva_agg
  sc_agg <- singlecell_result$gsva_agg

  common_pathways <- intersect(rownames(pb_agg), rownames(sc_agg))

  comparison <- data.frame(
    pathway = common_pathways,
    pb_HSD17B7_neg = pb_agg[common_pathways, "HSD17B7-"],
    pb_HSD17B7_pos = pb_agg[common_pathways, "HSD17B7+"],
    pb_diff = pb_agg[common_pathways, "HSD17B7+"] - pb_agg[common_pathways, "HSD17B7-"],
    sc_HSD17B7_neg = sc_agg[common_pathways, "HSD17B7-"],
    sc_HSD17B7_pos = sc_agg[common_pathways, "HSD17B7+"],
    sc_diff = sc_agg[common_pathways, "HSD17B7+"] - sc_agg[common_pathways, "HSD17B7-"],
    stringsAsFactors = FALSE
  )
  comparison$direction_match <- sign(comparison$pb_diff) == sign(comparison$sc_diff)

  write.csv(comparison, file.path(output_dir, "gsva_comparison_pb_vs_sc.csv"), row.names = FALSE)

  cat("\n  Direction comparison (HSD17B7+ vs HSD17B7-):\n")
  for (i in seq_len(nrow(comparison))) {
    dir_pb <- ifelse(comparison$pb_diff[i] > 0, "UP", "DOWN")
    dir_sc <- ifelse(comparison$sc_diff[i] > 0, "UP", "DOWN")
    match_str <- ifelse(comparison$direction_match[i], "MATCH", "MISMATCH")
    cat(sprintf("  %-45s PB: %s  SC: %s  [%s]\n",
                comparison$pathway[i], dir_pb, dir_sc, match_str))
  }

  n_match <- sum(comparison$direction_match)
  cat(sprintf("\n  %d/%d pathways agree in direction (%.0f%%)\n",
              n_match, nrow(comparison), 100 * n_match / nrow(comparison)))

  # Side-by-side heatmap
  cat("\n  Generating side-by-side comparison heatmap...\n")
  combined <- cbind(
    pb_agg[common_pathways, ],
    sc_agg[common_pathways, ]
  )
  colnames(combined) <- c("PB: HSD17B7-", "PB: HSD17B7+",
                           "SC: HSD17B7-", "SC: HSD17B7+")

  data_range <- max(abs(combined), na.rm = TRUE)
  if (data_range < 0.01) data_range <- 0.5

  col_fun <- colorRamp2(
    c(-data_range, 0, data_range),
    c("#313695", "white", "#d73027")
  )

  ht_cmp <- Heatmap(
    combined,
    name = "GSVA\nscore",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    column_split = factor(c("Pseudo-bulk", "Pseudo-bulk", "Single-cell", "Single-cell"),
                          levels = c("Pseudo-bulk", "Single-cell")),
    rect_gp = gpar(col = "grey90", lwd = 0.5),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", combined[i, j]), x, y, gp = gpar(fontsize = 7))
    }
  )

  png(file.path(figures_dir, "gsva_comparison_heatmap.png"), width = 14*300, height = 8*300, res = 300)
  draw(ht_cmp, padding = unit(c(2, 2, 2, 2), "cm"))
  dev.off()

  pdf(file.path(output_dir, "gsva_comparison_heatmap.pdf"), width = 14, height = 8)
  draw(ht_cmp, padding = unit(c(2, 2, 2, 2), "cm"))
  dev.off()
}

# Also copy the primary (pseudo-bulk) heatmap as the default for backwards compat
if (!is.null(pseudobulk_result)) {
  file.copy(file.path(figures_dir, "gsva_heatmap_pseudobulk.png"),
            file.path(figures_dir, "gsva_heatmap.png"), overwrite = TRUE)
  file.copy(file.path(output_dir, "gsva_heatmap_pseudobulk.pdf"),
            file.path(output_dir, "gsva_heatmap.pdf"), overwrite = TRUE)
}

cat("\n=== scRNA GSVA complete ===\n")
cat("Outputs:\n")
if (!is.null(pseudobulk_result)) {
  cat("  Pseudo-bulk: figures/gsva_heatmap_pseudobulk.png\n")
}
if (!is.null(singlecell_result)) {
  cat("  Single-cell: figures/gsva_heatmap_singlecell.png\n")
}
if (!is.null(pseudobulk_result) && !is.null(singlecell_result)) {
  cat("  Comparison:  figures/gsva_comparison_heatmap.png\n")
  cat("  CSV:         outputs/gsva_comparison_pb_vs_sc.csv\n")
}
