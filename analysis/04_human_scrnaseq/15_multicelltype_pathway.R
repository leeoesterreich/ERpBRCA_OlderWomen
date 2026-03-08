#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/15_multicelltype_pathway.R
# HALLMARK/BIOCARTA pathway enrichment per cell type
# Generates Figure 7 Panels B/C - dual pathway activity heatmaps
#
# Usage:
#   Rscript 15_multicelltype_pathway.R [--mode=curated|divergent] [--n-pathways=12]
#
# Modes:
#   curated   - Use hardcoded curated pathway list from manuscript (default)
#   divergent - Select top N pathways by Young vs Elderly z-score difference
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/macrophage_seurat.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/multicelltype_pathway_scores.csv
#   - analysis/04_human_scrnaseq/figures/fig7bc_pathway_heatmaps.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(msigdbr)
  library(GSVA)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tibble)
})

# -----------------------------------------------------------------------------
# Parse command line arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
pathway_mode <- "curated"  # default
n_pathways <- 12  # default number of pathways per category in divergent mode

for (arg in args) {
  if (grepl("^--mode=", arg)) {
    pathway_mode <- sub("^--mode=", "", arg)
  }
  if (grepl("^--n-pathways=", arg)) {
    n_pathways <- as.integer(sub("^--n-pathways=", "", arg))
  }
}

if (!pathway_mode %in% c("curated", "divergent")) {
  stop("Invalid mode. Use --mode=curated or --mode=divergent")
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
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Multi-Cell-Type Pathway Enrichment (Figure 7C) ===\n")
cat(sprintf("Mode: %s\n", pathway_mode))
if (pathway_mode == "divergent") {
  cat(sprintf("Selecting top %d pathways per category by age divergence\n", n_pathways))
}

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_file <- file.path(output_dir, "macrophage_seurat.rds")
seurat_obj <- readRDS(seurat_file)

# Filter to Elderly and Young
seurat_ey <- subset(seurat_obj, subset = AgeGroup %in% c("Elderly", "Young"))
cat("  Cells:", ncol(seurat_ey), "\n")

# -----------------------------------------------------------------------------
# Step 2: Load HALLMARK and BIOCARTA gene sets
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading gene sets...\n")

# HALLMARK - use ALL pathways (50 total) for better z-score range
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
cat("  HALLMARK pathways loaded:", length(hallmark_list), "\n")

# BIOCARTA - use ALL pathways
biocarta_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
biocarta_list <- split(biocarta_sets$gene_symbol, biocarta_sets$gs_name)
cat("  BIOCARTA pathways loaded:", length(biocarta_list), "\n")

all_pathways <- c(hallmark_list, biocarta_list)
cat("  Total pathways:", length(all_pathways), "\n")

# -----------------------------------------------------------------------------
# Step 3: Create pseudo-bulk per cell type per age group
# -----------------------------------------------------------------------------
cat("\nStep 3: Creating pseudo-bulk profiles...\n")

# FIX: Use SCTransform data if available (matches manuscript methods)
# SCT provides better normalization for pseudo-bulk GSVA than CPM+log2
# Magnitude analysis showed CPM inflates estrogen pathway effects by 25x vs SCT
if ("SCT" %in% Assays(seurat_ey)) {
  cat("  Using SCTransform data (matches manuscript methods)\n")
  DefaultAssay(seurat_ey) <- "SCT"
  # SCT 'data' slot contains corrected log-normalized counts
  expr_data <- GetAssayData(seurat_ey, layer = "data")
} else {
  cat("  WARNING: SCT assay not found, falling back to CPM+log2\n")
  cat("  This may produce different effect magnitudes than the manuscript\n")
  DefaultAssay(seurat_ey) <- "RNA"
  expr_data <- GetAssayData(seurat_ey, layer = "counts")
}

# Create grouping variable
seurat_ey$ct_age <- paste0(seurat_ey$CellTypeAnnotSH, "_", seurat_ey$AgeGroup)
groups <- unique(seurat_ey$ct_age)

# Aggregate (sum for counts, mean for SCT normalized data)
use_sct <- DefaultAssay(seurat_ey) == "SCT"

pseudobulk <- sapply(groups, function(grp) {
  cells <- colnames(seurat_ey)[seurat_ey$ct_age == grp]
  if (length(cells) > 10) {
    if (use_sct) {
      # For SCT data, take mean (already normalized)
      Matrix::rowMeans(expr_data[, cells, drop = FALSE])
    } else {
      # For counts, sum then normalize
      Matrix::rowSums(expr_data[, cells, drop = FALSE])
    }
  } else {
    rep(NA, nrow(expr_data))
  }
})

# Remove NA columns
pseudobulk <- pseudobulk[, !apply(pseudobulk, 2, function(x) all(is.na(x)))]

# Normalize only if using counts (SCT is already normalized)
if (use_sct) {
  pseudobulk_log <- pseudobulk  # SCT data is already log-normalized
  cat("  SCT data: no additional normalization needed\n")
} else {
  pseudobulk_cpm <- sweep(pseudobulk, 2, colSums(pseudobulk, na.rm = TRUE), "/") * 1e6
  pseudobulk_log <- log2(pseudobulk_cpm + 1)
  cat("  Applied CPM + log2 normalization\n")
}

cat("  Pseudo-bulk matrix:", nrow(pseudobulk_log), "genes x", ncol(pseudobulk_log), "groups\n")

# -----------------------------------------------------------------------------
# Step 4: Run GSVA
# -----------------------------------------------------------------------------
cat("\nStep 4: Running GSVA...\n")

gsva_result <- gsva(
  gsvaParam(
    as.matrix(pseudobulk_log),
    all_pathways,
    kcdf = "Gaussian",
    maxDiff = TRUE
  )
)

cat("  GSVA result:", nrow(gsva_result), "pathways x", ncol(gsva_result), "groups\n")

# Per-pathway z-score normalization (row-wise) to expand range
cat("  Applying per-pathway z-score normalization...\n")
gsva_result <- t(scale(t(gsva_result)))
cat("  Z-score range:", round(min(gsva_result, na.rm = TRUE), 2), "to",
    round(max(gsva_result, na.rm = TRUE), 2), "\n")

# Save normalized scores
gsva_df <- as.data.frame(gsva_result) %>%
  tibble::rownames_to_column("pathway")
fwrite(gsva_df, file.path(output_dir, "multicelltype_pathway_scores.csv"))
cat("  Saved multicelltype_pathway_scores.csv\n")

# -----------------------------------------------------------------------------
# Step 5: Generate Figure 7 B/C - Dual Heatmaps (HALLMARK + BIOCARTA)
# -----------------------------------------------------------------------------
cat("\nStep 5: Generating Figure 7 B/C dual heatmaps...\n")

# Transpose: rows = cell_type_age, columns = pathways
gsva_t <- t(gsva_result)

# Create row annotations for cell type categories
get_category <- function(ct_age) {
  ct <- gsub("_(Elderly|Young)$", "", ct_age)
  case_when(
    grepl("Tcells|NK|Bcells|Plasma", ct) ~ "Lymphocyte",
    grepl("Macro|Mono|DC|Myeloid", ct) ~ "Myeloid",
    grepl("Epithelial|Cancer|Luminal|Basal", ct) ~ "Epithelial",
    grepl("CAF|PVL|Endo|Fibro", ct) ~ "Stromal",
    TRUE ~ "Other"
  )
}

row_categories <- data.frame(
  Category = sapply(rownames(gsva_t), get_category),
  row.names = rownames(gsva_t)
)

# Define row order: group by category, then by cell type (Younger before Older)
category_order <- c("Lymphocyte", "Myeloid", "Epithelial", "Stromal", "Other")
row_order <- rownames(gsva_t)[order(
  match(row_categories$Category, category_order),
  gsub("_(Elderly|Young)$", "", rownames(gsva_t)),
  grepl("_Elderly$", rownames(gsva_t))
)]
gsva_t <- gsva_t[row_order, ]
row_categories <- row_categories[row_order, , drop = FALSE]

# Rename row labels: "Young" -> "Younger", "Elderly" -> "Older" (match manuscript)
rownames(gsva_t) <- gsub("_Young$", "_Younger", rownames(gsva_t))
rownames(gsva_t) <- gsub("_Elderly$", "_Older", rownames(gsva_t))
rownames(row_categories) <- rownames(gsva_t)

# Calculate gaps for category separation
category_counts <- table(row_categories$Category)[category_order]
category_counts <- category_counts[!is.na(category_counts) & category_counts > 0]
gaps_row <- cumsum(category_counts)[-length(category_counts)]

# Split into HALLMARK and BIOCARTA
hallmark_cols <- grep("^HALLMARK_", colnames(gsva_t), value = TRUE)
biocarta_cols <- grep("^BIOCARTA_", colnames(gsva_t), value = TRUE)

# -----------------------------------------------------------------------------
# Pathway selection based on mode
# -----------------------------------------------------------------------------
if (pathway_mode == "divergent") {
  cat("\nStep 5a: Selecting divergent pathways (Younger vs Older)...\n")

  # Function to select top N pathways by age divergence
  select_divergent_pathways <- function(mat, pathway_cols, n_select) {
    # Get Younger and Older row indices
    younger_rows <- grep("_Younger$", rownames(mat))
    older_rows <- grep("_Older$", rownames(mat))

    if (length(younger_rows) == 0 || length(older_rows) == 0) {
      warning("Cannot find Younger/Older rows for divergence calculation")
      return(pathway_cols[1:min(n_select, length(pathway_cols))])
    }

    # Calculate mean z-score per pathway for each age group
    divergence <- sapply(pathway_cols, function(pw) {
      younger_mean <- mean(mat[younger_rows, pw], na.rm = TRUE)
      older_mean <- mean(mat[older_rows, pw], na.rm = TRUE)
      abs(younger_mean - older_mean)  # Absolute difference
    })

    # Sort by divergence and take top N
    sorted_pws <- names(sort(divergence, decreasing = TRUE))
    selected <- head(sorted_pws, n_select)

    cat(sprintf("    Top divergent pathways (showing top 5 of %d):\n", length(selected)))
    for (i in 1:min(5, length(selected))) {
      pw <- selected[i]
      cat(sprintf("      %d. %s (divergence: %.2f)\n", i, pw, divergence[pw]))
    }

    return(selected)
  }

  hallmark_keep <- select_divergent_pathways(gsva_t, hallmark_cols, n_pathways)
  biocarta_keep <- select_divergent_pathways(gsva_t, biocarta_cols, n_pathways)

  cat(sprintf("  Selected %d HALLMARK and %d BIOCARTA divergent pathways\n",
              length(hallmark_keep), length(biocarta_keep)))

} else {
  # CURATED mode (default) - manuscript shows ~12 pathways per panel
  # These are the pathways visible in the manuscript Figure 7B/C
  # FIX: Added TGF_BETA_SIGNALING - manuscript explicitly claims "TGFβ signaling enriched in older"
  hallmark_curated <- c(
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_HYPOXIA",
    "HALLMARK_APOPTOSIS"
  )

  biocarta_curated <- c(
    "BIOCARTA_INFLAM_PATHWAY",
    "BIOCARTA_IL6_PATHWAY",
    "BIOCARTA_IL2_PATHWAY",
    "BIOCARTA_NFKB_PATHWAY",
    "BIOCARTA_TNFR1_PATHWAY",
    "BIOCARTA_DEATH_PATHWAY",
    "BIOCARTA_FAS_PATHWAY",
    "BIOCARTA_CASPASE_PATHWAY",
    "BIOCARTA_P53_PATHWAY",
    "BIOCARTA_CELLCYCLE_PATHWAY",
    "BIOCARTA_G1_PATHWAY",
    "BIOCARTA_G2_PATHWAY"
  )

  # Filter to curated pathways (keep order)
  hallmark_keep <- intersect(hallmark_curated, hallmark_cols)
  biocarta_keep <- intersect(biocarta_curated, biocarta_cols)

  # FIX: Warn about missing curated pathways (prevents silent dropping)
  hallmark_missing <- setdiff(hallmark_curated, hallmark_cols)
  biocarta_missing <- setdiff(biocarta_curated, biocarta_cols)

  if (length(hallmark_missing) > 0) {
    warning(sprintf("Missing HALLMARK pathways (not in GSVA output): %s",
                    paste(hallmark_missing, collapse = ", ")))
  }
  if (length(biocarta_missing) > 0) {
    warning(sprintf("Missing BIOCARTA pathways (not in GSVA output): %s",
                    paste(biocarta_missing, collapse = ", ")))
  }

  cat("  Using", length(hallmark_keep), "HALLMARK and", length(biocarta_keep), "BIOCARTA curated pathways\n")
}

# Assert that we have at least some pathways
if (length(hallmark_keep) == 0) {
  stop("ERROR: No HALLMARK pathways matched. Check pathway names or GSVA output.")
}
if (length(biocarta_keep) == 0) {
  stop("ERROR: No BIOCARTA pathways matched. Check pathway names or GSVA output.")
}

gsva_hallmark <- gsva_t[, hallmark_keep, drop = FALSE]
gsva_biocarta <- gsva_t[, biocarta_keep, drop = FALSE]

# Clean column names for display
colnames(gsva_hallmark) <- gsub("^HALLMARK_", "", colnames(gsva_hallmark))
colnames(gsva_hallmark) <- gsub("_", " ", colnames(gsva_hallmark))
colnames(gsva_biocarta) <- gsub("^BIOCARTA_", "", colnames(gsva_biocarta))
colnames(gsva_biocarta) <- gsub("_", " ", colnames(gsva_biocarta))

# Annotation colors (manuscript-matching saturated colors)
ann_colors <- list(
  Category = c(
    Lymphocyte = "#4DAF4A",
    Myeloid = "#E41A1C",
    Epithelial = "#377EB8",
    Stromal = "#984EA3",
    Other = "grey70"
  )
)

# Manuscript color palette: blue-cyan -> near-white -> orange-red
# Sampled from manuscript: #8BC5E4, #BCDEEE -> #EDF3F1 -> #F6BD6A, #E75321
hallmark_col_fun <- colorRamp2(
  c(-3, -1.5, 0, 1.5, 3),
  c("#3E5CA8", "#8BC5E4", "#EDF3F1", "#F6BD6A", "#E75321")
)
biocarta_col_fun <- colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("#3E5CA8", "#8BC5E4", "#EDF3F1", "#F6BD6A", "#E75321")
)

# Row split by category (manuscript-style left labels)
row_split <- factor(
  row_categories$Category,
  levels = category_order
)

# Left annotation bar
left_anno <- rowAnnotation(
  Category = row_categories$Category,
  col = ann_colors,
  show_annotation_name = FALSE,
  width = unit(3, "mm"),
  border = FALSE
)

# Helper function for consistent heatmap styling
make_manuscript_heatmap <- function(mat, title_text, col_fun, range_vals) {
  Heatmap(
    mat,
    name = "z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    row_split = row_split,
    row_gap = unit(2, "mm"),
    row_title_side = "left",
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 10, fontface = "plain"),
    left_annotation = left_anno,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 9),
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 10),  # Larger font for readable labels
    column_title = title_text,
    column_title_gp = gpar(fontsize = 16, fontface = "plain"),
    column_dend_height = unit(10, "mm"),
    border = FALSE,
    rect_gp = gpar(col = "#3B3B3B", lwd = 0.5),  # Thin dark cell borders
    width = unit(55, "mm"),   # Adjusted for ~12 columns
    height = unit(110, "mm"),
    heatmap_legend_param = list(
      direction = "horizontal",
      title = "Pathway activity\n(z-score)",
      at = c(range_vals[1], 0, range_vals[2]),
      labels = c(as.character(range_vals[1]), "0", as.character(range_vals[2])),
      legend_width = unit(30, "mm"),
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 9)
    )
  )
}

# Create heatmaps with manuscript styling
# Title case "Hallmark" not "HALLMARK"
ht_hallmark <- make_manuscript_heatmap(gsva_hallmark, "Hallmark", hallmark_col_fun, c(-3, 3))
ht_biocarta <- make_manuscript_heatmap(gsva_biocarta, "BIOCARTA", biocarta_col_fun, c(-4, 4))

# Determine output filename based on mode
mode_suffix <- ifelse(pathway_mode == "divergent", "_divergent", "")

# Draw side-by-side with legend at bottom (manuscript layout)
png_file <- file.path(figures_dir, sprintf("fig7bc_pathway_heatmaps%s.png", mode_suffix))
png(png_file, width = 3000, height = 2400, res = 300, bg = "white")
draw(
  ht_hallmark + ht_biocarta,
  heatmap_legend_side = "bottom",
  merge_legends = TRUE,
  padding = unit(c(10, 15, 10, 10), "mm")
)
dev.off()
cat(sprintf("  Saved %s\n", basename(png_file)))

# SVG for vector graphics
svg_file <- file.path(figures_dir, sprintf("fig7bc_pathway_heatmaps%s.svg", mode_suffix))
svg(svg_file, width = 10, height = 8)
draw(
  ht_hallmark + ht_biocarta,
  heatmap_legend_side = "bottom",
  merge_legends = TRUE,
  padding = unit(c(10, 15, 10, 10), "mm")
)
dev.off()
cat(sprintf("  Saved %s\n", basename(svg_file)))

cat(sprintf("\n=== Multi-cell-type pathway analysis complete (mode: %s) ===\n", pathway_mode))
