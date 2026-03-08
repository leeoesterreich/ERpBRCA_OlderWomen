#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/06_run_gsva.R
# GSVA on single-cell data using pseudo-bulk approach
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/gsva_pseudobulk.rds
#   - analysis/04_human_scrnaseq/outputs/gsva_heatmap.pdf
#   - analysis/04_human_scrnaseq/figures/gsva_heatmap.png

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
  library(pheatmap)  # For supplementary per-sample heatmap
  library(tibble)
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
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== scRNA-seq GSVA (Pseudo-bulk) ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))
cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Create pseudo-bulk per patient
# -----------------------------------------------------------------------------
cat("\nStep 2: Creating pseudo-bulk profiles...\n")

# Aggregate counts per patient
DefaultAssay(seurat_obj) <- "RNA"

# Get raw counts
counts <- GetAssayData(seurat_obj, slot = "counts")

# Aggregate by patient (orig.ident)
patients <- unique(seurat_obj$orig.ident)
pseudobulk <- sapply(patients, function(pt) {
  cells <- colnames(seurat_obj)[seurat_obj$orig.ident == pt]
  rowSums(counts[, cells, drop = FALSE])
})

# Normalize (CPM + log)
pseudobulk_cpm <- sweep(pseudobulk, 2, colSums(pseudobulk), "/") * 1e6
pseudobulk_log <- log2(pseudobulk_cpm + 1)

cat("  Pseudo-bulk matrix:", nrow(pseudobulk_log), "genes x", ncol(pseudobulk_log), "samples\n")

# -----------------------------------------------------------------------------
# Step 3: Load gene sets (manuscript Figure 4E - 11 pathways in fixed order)
# -----------------------------------------------------------------------------
cat("\nStep 3: Loading gene sets (manuscript order)...\n")

# Get all gene sets
all_sets <- msigdbr(species = "Homo sapiens")
all_list <- split(all_sets$gene_symbol, all_sets$gs_name)

# Helper function to find pathway by pattern with verbose debugging
find_pathway <- function(pattern, verbose = TRUE) {
  hit <- grep(pattern, names(all_list), ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) {
    if (verbose) cat("  WARNING: No match for pattern '", pattern, "'\n", sep = "")
    return(NULL)
  }
  if (verbose) cat("  Found: ", hit[1], " (", length(all_list[[hit[1]]]), " genes)\n", sep = "")
  all_list[[hit[1]]]
}

# Manuscript Figure 4E: 11 pathways in display order
# Map display names to msigdbr gene set names
manuscript_display_order <- c(
  "Estrogen metabolism WP697",
  "Estrogen metabolism WP5276",
  "Estrogen receptor pathway",
  "Response to estrogen",
  "Estrogen signaling",
  "Estrogen receptor signaling pathway",
  "Estrogen dependent gene expression",
  "Positive regulation of ER signaling",
  "Estrogen response late",
  "Estrogen response early",
  "Cellular response to estrogen stimulus"
)

# Build pathway list with fallback patterns for robustness
cat("  Looking up pathways...\n")
estrogen_pathways <- list(
  "Estrogen metabolism WP697" = find_pathway("WP_ESTROGEN_METABOLISM"),
  "Estrogen metabolism WP5276" = find_pathway("WP_ESTROGEN_BIOSYNTHESIS"),
  "Estrogen receptor pathway" = {
    # BIOCARTA may have different naming conventions
    p <- find_pathway("BIOCARTA.*ESTROGEN", verbose = FALSE)
    if (is.null(p)) p <- find_pathway("ESTROGEN.*BIOCARTA", verbose = FALSE)
    if (is.null(p)) p <- find_pathway("PID_ER_PATHWAY", verbose = FALSE)  # Alternative
    if (is.null(p)) cat("  WARNING: Estrogen receptor pathway not found\n")
    p
  },
  "Response to estrogen" = find_pathway("GOBP_RESPONSE_TO_ESTROGEN$"),
  "Estrogen signaling" = {
    p <- find_pathway("WP_ESTROGEN_SIGNALING", verbose = FALSE)
    if (is.null(p)) p <- find_pathway("KEGG_ESTROGEN_SIGNALING", verbose = FALSE)
    if (is.null(p)) cat("  WARNING: Estrogen signaling pathway not found\n")
    p
  },
  "Estrogen receptor signaling pathway" = find_pathway("GOBP_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY"),
  "Estrogen dependent gene expression" = find_pathway("REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION"),
  "Positive regulation of ER signaling" = find_pathway("GOBP_POSITIVE_REGULATION_OF_INTRACELLULAR_ESTROGEN"),
  "Estrogen response late" = find_pathway("HALLMARK_ESTROGEN_RESPONSE_LATE"),
  "Estrogen response early" = find_pathway("HALLMARK_ESTROGEN_RESPONSE_EARLY"),
  "Cellular response to estrogen stimulus" = find_pathway("GOBP_CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS")
)

# Remove any NULL entries
estrogen_pathways <- estrogen_pathways[!sapply(estrogen_pathways, is.null)]

# Add LI_ESTROGENE E2 response signatures from local GMT files
gmt_dir <- file.path(project_root, "data", "gmt")
parse_gmt <- function(path) {
  line <- readLines(path, n = 1)
  fields <- strsplit(line, "\t")[[1]]
  fields[-(1:2)]  # skip name and URL
}
estrogen_pathways[["LI EstroGene early E2 response up"]] <- parse_gmt(file.path(gmt_dir, "LI_ESTROGENE_EARLY_E2_RESPONSE_UP.v2025.1.Hs.gmt"))
estrogen_pathways[["LI EstroGene late E2 response up"]] <- parse_gmt(file.path(gmt_dir, "LI_ESTROGENE_LATE_E2_RESPONSE_UP.v2025.1.Hs.gmt"))

cat("  Gene sets:", length(estrogen_pathways), "\n")

# -----------------------------------------------------------------------------
# Step 4: Run GSVA
# -----------------------------------------------------------------------------
cat("\nStep 4: Running GSVA...\n")

gsva_result <- gsva(
  gsvaParam(
    as.matrix(pseudobulk_log),
    estrogen_pathways,
    kcdf = "Gaussian",
    maxDiff = TRUE
  )
)

cat("  Result:", nrow(gsva_result), "pathways x", ncol(gsva_result), "samples\n")

# -----------------------------------------------------------------------------
# Step 5: Aggregate by HSD17B7 status and visualize (match manuscript Fig 4E)
# -----------------------------------------------------------------------------
cat("\nStep 5: Aggregating by HSD17B7 status and visualizing...\n")

# Get HSD17B7 expression per sample
hsd17b7_expr <- pseudobulk_log["HSD17B7", colnames(gsva_result)]

# Split samples by median HSD17B7 expression
hsd17b7_median <- median(hsd17b7_expr, na.rm = TRUE)
hsd17b7_group <- ifelse(hsd17b7_expr > hsd17b7_median, "HSD17B7+", "HSD17B7-")
cat("  HSD17B7 median:", round(hsd17b7_median, 2), "\n")
cat("  HSD17B7- samples:", sum(hsd17b7_group == "HSD17B7-"), "\n")
cat("  HSD17B7+ samples:", sum(hsd17b7_group == "HSD17B7+"), "\n")

# Aggregate GSVA scores by HSD17B7 status (mean per group)
gsva_aggregated <- sapply(c("HSD17B7-", "HSD17B7+"), function(g) {
  cols <- names(hsd17b7_group)[hsd17b7_group == g]
  rowMeans(gsva_result[, cols, drop = FALSE], na.rm = TRUE)
})

# Use RAW GSVA aggregated scores - DO NOT scale per-row
# The manuscript shows natural gradient where HSD17B7- is near 0 (white)
# and HSD17B7+ shows actual pathway activation (red) or suppression (blue)
# Row-scaling forces binary ±4 extremes which looks wrong
gsva_scaled <- gsva_aggregated

# Clamp to -4/+4 range for visualization (but keep natural variation)
gsva_scaled <- pmax(pmin(gsva_scaled, 4), -4)

# Ensure manuscript row order (only keep rows that exist)
valid_rows <- intersect(names(estrogen_pathways), rownames(gsva_scaled))
gsva_scaled <- gsva_scaled[valid_rows, , drop = FALSE]

cat("  Aggregated matrix:", nrow(gsva_scaled), "pathways x", ncol(gsva_scaled), "groups\n")
cat("  Data range:", round(min(gsva_scaled, na.rm = TRUE), 2), "to",
    round(max(gsva_scaled, na.rm = TRUE), 2), "\n")

# -----------------------------------------------------------------------------
# Create manuscript-matching heatmap using ComplexHeatmap
# Manuscript style: row labels LEFT, column labels 45°, no cell borders,
# compact cells, blue-white-red colorbar labeled "Pathway activity"
# -----------------------------------------------------------------------------

# Color scale - use actual data range, not fixed -4/+4
# Raw GSVA scores are typically in range -0.5 to +0.5
data_range <- max(abs(gsva_scaled), na.rm = TRUE)
cat("  Color scale range: -", round(data_range, 2), " to +", round(data_range, 2), "\n", sep = "")

col_fun <- colorRamp2(
  c(-data_range, -data_range/2, 0, data_range/2, data_range),
  c("#313695", "#74add1", "white", "#fdae61", "#d73027")
)

# Create heatmap with manuscript styling
ht <- Heatmap(
  gsva_scaled,
  name = "Pathway\nactivity",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",              # Row labels on LEFT (manuscript style)
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,                # Column labels rotated 45°
  column_names_gp = gpar(fontsize = 12),
  column_names_side = "bottom",
  rect_gp = gpar(col = NA),             # NO cell borders
  width = unit(2, "cm"),                # Compact cells
  height = unit(nrow(gsva_scaled) * 0.5, "cm"),
  heatmap_legend_param = list(
    title = "Pathway\nactivity",
    at = c(-round(data_range, 1), 0, round(data_range, 1)),
    labels = c(as.character(-round(data_range, 1)), "0", as.character(round(data_range, 1))),
    legend_height = unit(3, "cm")
  )
)

# PDF output
pdf(file.path(output_dir, "gsva_heatmap.pdf"), width = 10, height = 6)
draw(ht, padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()

# PNG for validation pipeline
png(file.path(figures_dir, "gsva_heatmap.png"), width = 10*300, height = 6*300, res = 300)
draw(ht, padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()

# SVG for vector graphics
svg(file.path(figures_dir, "gsva_heatmap.svg"), width = 10, height = 6)
draw(ht, padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()

# Also save per-sample version for supplementary
cat("\n  Saving per-sample heatmap (supplementary)...\n")
age_groups <- seurat_obj@meta.data %>%
  select(orig.ident, AgeGroup) %>%
  distinct()
rownames(age_groups) <- age_groups$orig.ident

ann_col <- data.frame(
  AgeGroup = age_groups[colnames(gsva_result), "AgeGroup"],
  HSD17B7 = hsd17b7_group[colnames(gsva_result)],
  row.names = colnames(gsva_result)
)

ann_colors <- list(
  AgeGroup = c(Young = "#4DAF4A", MidAge = "#377EB8", Elderly = "#E41A1C"),
  HSD17B7 = c("HSD17B7-" = "#74add1", "HSD17B7+" = "#d73027")
)

png(file.path(figures_dir, "gsva_heatmap_per_sample.png"), width = 12*300, height = 8*300, res = 300)
pheatmap(
  gsva_result,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  scale = "row",
  show_colnames = TRUE,
  main = NA
)
dev.off()

# -----------------------------------------------------------------------------
# Step 6: Save
# -----------------------------------------------------------------------------
cat("\nStep 6: Saving...\n")

saveRDS(gsva_result, file.path(output_dir, "gsva_pseudobulk.rds"))
saveRDS(gsva_scaled, file.path(output_dir, "gsva_aggregated_by_hsd17b7.rds"))
saveRDS(pseudobulk_log, file.path(output_dir, "pseudobulk_log2cpm.rds"))

# Save aggregated scores as CSV for easy inspection
write.csv(gsva_scaled, file.path(output_dir, "gsva_aggregated_by_hsd17b7.csv"))

cat("\n=== scRNA GSVA complete ===\n")
cat("Main figure: figures/gsva_heatmap.png (manuscript Fig 4E format)\n")
cat("Supplementary: figures/gsva_heatmap_per_sample.png (per-sample view)\n")
