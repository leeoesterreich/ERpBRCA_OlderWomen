#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/15_multicelltype_pathway.R
# HALLMARK/BIOCARTA pathway enrichment per cell type with per-patient pseudo-bulk
# Generates Figure 7 Panels B/C - dual pathway activity heatmaps
#
# Approach: For each cell type, create pseudo-bulk per PATIENT (not per age group),
# run GSVA across all patients within that cell type, then compare Young vs Elderly
# with Welch t-test tests. This preserves patient-level replication.
#
# Usage:
#   Rscript 15_multicelltype_pathway.R [--mode=curated|divergent] [--n-pathways=12]
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/macrophage_seurat.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/multicelltype_pathway_scores.csv
#   - analysis/04_human_scrnaseq/outputs/multicelltype_pathway_stats.csv
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
pathway_mode <- "curated"
n_pathways <- 12

for (arg in args) {
  if (grepl("^--mode=", arg)) pathway_mode <- sub("^--mode=", "", arg)
  if (grepl("^--n-pathways=", arg)) n_pathways <- as.integer(sub("^--n-pathways=", "", arg))
}

if (!pathway_mode %in% c("curated", "divergent")) {
  stop("Invalid mode. Use --mode=curated or --mode=divergent")
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("--file=", "", file_arg))))
  return(getwd())
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "outputs")
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Multi-Cell-Type Pathway Enrichment (Figure 7B/C) ===\n")
cat(sprintf("Mode: %s\n", pathway_mode))
cat("Approach: per-patient pseudo-bulk GSVA with Welch's t-test (var.equal=FALSE)\n")

# -----------------------------------------------------------------------------
# Step 1: Load data — use ALL age groups (Young + MidAge + Elderly)
# -----------------------------------------------------------------------------
cat("\nStep 1: Loading data...\n")

seurat_file <- file.path(output_dir, "seurat_annotated.rds")
seurat_obj <- readRDS(seurat_file)

# Fix tab-embedded gene names from Xu atlas (e.g., "ESR1\tESR1" → "ESR1")
gene_names <- rownames(seurat_obj)
if (any(grepl("\t", gene_names))) {
  cat("  Cleaning tab-embedded gene names...\n")
  clean_names <- sapply(strsplit(gene_names, "\t"), function(x) x[length(x)])
  clean_names <- make.unique(clean_names)
  for (assay_name in Assays(seurat_obj)) {
    assay_obj <- seurat_obj[[assay_name]]
    rownames(assay_obj) <- clean_names[match(rownames(assay_obj), gene_names)]
    seurat_obj[[assay_name]] <- assay_obj
  }
  cat("  Cleaned", sum(grepl("\t", gene_names)), "gene names\n")
}

# Keep all age groups for GSVA (more samples = better z-score estimation)
# Statistical tests will compare Young vs Elderly
cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  Age groups:\n")
print(table(seurat_obj$AgeGroup))
cat("  Cell types:\n")
print(table(seurat_obj$CellTypeAnnotSH))

# -----------------------------------------------------------------------------
# Step 2: Load gene sets
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading gene sets...\n")

hallmark_sets <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H")
)
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
cat("  HALLMARK pathways loaded:", length(hallmark_list), "\n")

biocarta_sets <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:BIOCARTA"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
)
biocarta_list <- split(biocarta_sets$gene_symbol, biocarta_sets$gs_name)
cat("  BIOCARTA pathways loaded:", length(biocarta_list), "\n")

gmt_dir <- file.path(normalizePath(file.path(script_dir, "../..")), "data", "gmt")
parse_gmt <- function(path) {
  line <- readLines(path, n = 1)
  fields <- strsplit(line, "\t")[[1]]
  fields[-(1:2)]
}
li_estrogene_list <- list(
  LI_ESTROGENE_EARLY_E2_RESPONSE_UP = parse_gmt(file.path(gmt_dir, "LI_ESTROGENE_EARLY_E2_RESPONSE_UP.v2025.1.Hs.gmt")),
  LI_ESTROGENE_LATE_E2_RESPONSE_UP = parse_gmt(file.path(gmt_dir, "LI_ESTROGENE_LATE_E2_RESPONSE_UP.v2025.1.Hs.gmt"))
)
cat("  LI_ESTROGENE pathways loaded:", length(li_estrogene_list), "\n")

all_pathways <- c(hallmark_list, biocarta_list, li_estrogene_list)
cat("  Total pathways:", length(all_pathways), "\n")

# -----------------------------------------------------------------------------
# Step 3: Per-patient pseudo-bulk within each cell type
# Run GSVA per cell type (patients as columns)
# -----------------------------------------------------------------------------
cat("\nStep 3: Creating per-patient pseudo-bulk and running GSVA per cell type...\n")

min_cells_per_patient <- 10  # minimum cells for a patient to be included in a cell type
min_patients_per_group <- 2  # minimum patients per age group for statistical testing

if ("SCT" %in% Assays(seurat_obj)) {
  cat("  Using SCTransform data\n")
  DefaultAssay(seurat_obj) <- "SCT"
  expr_data <- GetAssayData(seurat_obj, layer = "data")
} else {
  cat("  WARNING: SCT not found, using RNA counts + CPM\n")
  DefaultAssay(seurat_obj) <- "RNA"
  expr_data <- GetAssayData(seurat_obj, layer = "counts")
}

# Clean tab-embedded gene names from expression matrix (Xu atlas artifact)
if (any(grepl("\t", rownames(expr_data)))) {
  cat("  Cleaning tab-embedded gene names from expression matrix...\n")
  clean_rn <- sapply(strsplit(rownames(expr_data), "\t"), function(x) x[length(x)])
  rownames(expr_data) <- make.unique(clean_rn)
  cat("  Expression matrix:", nrow(expr_data), "genes x", ncol(expr_data), "cells\n")
}

cell_types <- sort(unique(as.character(seurat_obj$CellTypeAnnotSH)))
patients <- sort(unique(as.character(seurat_obj$orig.ident)))

# Build patient -> age group mapping
patient_age <- seurat_obj@meta.data %>%
  select(orig.ident, AgeGroup) %>%
  distinct() %>%
  arrange(orig.ident)
patient_age_map <- setNames(as.character(patient_age$AgeGroup), as.character(patient_age$orig.ident))

# Store all results
all_gsva_scores <- list()  # per cell type: pathway x patient matrix
all_stats <- list()        # statistical test results

for (ct in cell_types) {
  cat(sprintf("\n  --- %s ---\n", ct))

  # Get cells for this cell type
  ct_cells <- colnames(seurat_obj)[seurat_obj$CellTypeAnnotSH == ct]
  ct_patients <- unique(as.character(seurat_obj$orig.ident[seurat_obj$CellTypeAnnotSH == ct]))

  # Create pseudo-bulk per patient
  pb_list <- list()
  for (pat in ct_patients) {
    pat_cells <- intersect(ct_cells, colnames(seurat_obj)[seurat_obj$orig.ident == pat])
    if (length(pat_cells) >= min_cells_per_patient) {
      pb_list[[pat]] <- Matrix::rowMeans(expr_data[, pat_cells, drop = FALSE])
    }
  }

  if (length(pb_list) < 3) {
    cat(sprintf("    Skipping: only %d patients with >= %d cells\n", length(pb_list), min_cells_per_patient))
    next
  }

  pb_mat <- do.call(cbind, pb_list)
  colnames(pb_mat) <- names(pb_list)

  # If using counts, normalize to CPM + log2
  if (DefaultAssay(seurat_obj) != "SCT") {
    pb_cpm <- sweep(pb_mat, 2, colSums(pb_mat), "/") * 1e6
    pb_mat <- log2(pb_cpm + 1)
  }

  cat(sprintf("    Pseudo-bulk: %d genes x %d patients\n", nrow(pb_mat), ncol(pb_mat)))

  # Run GSVA on this cell type's patient profiles
  gsva_ct <- tryCatch({
    gsva(gsvaParam(as.matrix(pb_mat), all_pathways, kcdf = "Gaussian", maxDiff = TRUE))
  }, error = function(e) {
    cat(sprintf("    WARNING: GSVA failed: %s\n", e$message))
    NULL
  })

  if (is.null(gsva_ct)) next

  all_gsva_scores[[ct]] <- gsva_ct

  # Statistical tests: Young vs Elderly for each pathway
  pat_ages <- patient_age_map[colnames(gsva_ct)]
  young_idx <- which(pat_ages == "Young")
  elderly_idx <- which(pat_ages == "Elderly")

  cat(sprintf("    Young patients: %d, Elderly patients: %d\n", length(young_idx), length(elderly_idx)))

  if (length(young_idx) >= min_patients_per_group && length(elderly_idx) >= min_patients_per_group) {
    for (pw in rownames(gsva_ct)) {
      young_scores <- gsva_ct[pw, young_idx]
      elderly_scores <- gsva_ct[pw, elderly_idx]

      # t-test (more powerful than Wilcoxon at small n)
      wt <- tryCatch(
        t.test(elderly_scores, young_scores, var.equal = FALSE),
        error = function(e) list(p.value = NA, statistic = NA)
      )

      all_stats[[length(all_stats) + 1]] <- data.frame(
        celltype = ct,
        pathway = pw,
        young_mean = mean(young_scores),
        elderly_mean = mean(elderly_scores),
        diff = mean(elderly_scores) - mean(young_scores),
        young_sd = sd(young_scores),
        elderly_sd = sd(elderly_scores),
        n_young = length(young_idx),
        n_elderly = length(elderly_idx),
        pvalue = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine stats and apply FDR correction
stats_df <- do.call(rbind, all_stats)

# FDR correction across ALL pathways (full correction)
stats_df$padj_all <- p.adjust(stats_df$pvalue, method = "BH")

# Also apply FDR correction ONLY to the curated pathways (reduced burden)
# This is the key test: ~24 pathways × ~13 cell types ≈ 312 tests vs 4000+
hallmark_curated <- c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS",
  "LI_ESTROGENE_EARLY_E2_RESPONSE_UP", "LI_ESTROGENE_LATE_E2_RESPONSE_UP"
)
biocarta_curated <- c(
  "BIOCARTA_INFLAM_PATHWAY", "BIOCARTA_IL6_PATHWAY", "BIOCARTA_IL2_PATHWAY",
  "BIOCARTA_NFKB_PATHWAY", "BIOCARTA_TNFR1_PATHWAY", "BIOCARTA_DEATH_PATHWAY",
  "BIOCARTA_FAS_PATHWAY", "BIOCARTA_CASPASE_PATHWAY", "BIOCARTA_P53_PATHWAY",
  "BIOCARTA_CELLCYCLE_PATHWAY", "BIOCARTA_G1_PATHWAY", "BIOCARTA_G2_PATHWAY"
)
curated_set <- c(hallmark_curated, biocarta_curated)
stats_df$is_curated <- stats_df$pathway %in% curated_set
curated_idx <- which(stats_df$is_curated)

# Per-cell-type BH correction (each cell type is a separate hypothesis space)
stats_df$padj <- NA_real_
for (ct in unique(stats_df$celltype)) {
  ct_curated <- which(stats_df$celltype == ct & stats_df$is_curated)
  if (length(ct_curated) > 0) {
    stats_df$padj[ct_curated] <- p.adjust(stats_df$pvalue[ct_curated], method = "BH")
  }
}
stats_df$significant <- !is.na(stats_df$padj) & stats_df$padj < 0.05

cat(sprintf("\n--- Statistical summary (t-test, curated FDR) ---\n"))
cat(sprintf("Total tests (all pathways): %d\n", nrow(stats_df)))
cat(sprintf("Curated pathway tests: %d\n", length(curated_idx)))
cat(sprintf("Significant curated (FDR < 0.05): %d\n", sum(stats_df$significant, na.rm = TRUE)))
cat(sprintf("Significant curated (FDR < 0.10): %d\n",
    sum(!is.na(stats_df$padj) & stats_df$padj < 0.10, na.rm = TRUE)))

# Also report nominal p < 0.05 (uncorrected) for curated
nominal_sig <- stats_df[curated_idx, ]
nominal_sig <- nominal_sig[!is.na(nominal_sig$pvalue) & nominal_sig$pvalue < 0.05, ]
cat(sprintf("Nominally significant curated (p < 0.05 uncorrected): %d\n", nrow(nominal_sig)))
if (nrow(nominal_sig) > 0) {
  nominal_sig <- nominal_sig[order(nominal_sig$pvalue), ]
  cat("\nNominally significant curated pathway x cell type (raw p < 0.05):\n")
  for (i in 1:min(20, nrow(nominal_sig))) {
    cat(sprintf("  %s | %s | diff=%+.3f | p=%.4f | padj=%.4f\n",
        nominal_sig$celltype[i], gsub("HALLMARK_|BIOCARTA_", "", nominal_sig$pathway[i]),
        nominal_sig$diff[i], nominal_sig$pvalue[i], nominal_sig$padj[i]))
  }
}

# Show FDR-significant results
sig <- stats_df[!is.na(stats_df$padj) & stats_df$padj < 0.10, ]
if (nrow(sig) > 0) {
  sig <- sig[order(sig$padj), ]
  cat("\nFDR-significant curated pathway x cell type (FDR < 0.10):\n")
  for (i in 1:min(20, nrow(sig))) {
    cat(sprintf("  %s | %s | diff=%+.3f | p=%.4f | padj=%.4f\n",
        sig$celltype[i], gsub("HALLMARK_|BIOCARTA_", "", sig$pathway[i]),
        sig$diff[i], sig$pvalue[i], sig$padj[i]))
  }
} else {
  cat("\nNo FDR-significant results even with curated-only correction.\n")
}

# Save stats
fwrite(stats_df, file.path(output_dir, "multicelltype_pathway_stats.csv"))
cat("  Saved multicelltype_pathway_stats.csv\n")

# -----------------------------------------------------------------------------
# Step 4: Build heatmap matrix — mean GSVA score per celltype x age group
# Now using properly computed per-patient scores, averaged per age group
# -----------------------------------------------------------------------------
cat("\nStep 4: Building heatmap matrix from per-patient scores...\n")

heatmap_rows <- list()
for (ct in names(all_gsva_scores)) {
  gsva_ct <- all_gsva_scores[[ct]]
  pat_ages <- patient_age_map[colnames(gsva_ct)]

  for (ag in c("Young", "Elderly")) {
    idx <- which(pat_ages == ag)
    if (length(idx) >= 1) {
      mean_scores <- rowMeans(gsva_ct[, idx, drop = FALSE])
      heatmap_rows[[paste0(ct, "_", ag)]] <- mean_scores
    }
  }
}

# Build matrix: rows = celltype_age, columns = pathways
gsva_mat <- do.call(rbind, heatmap_rows)

# Z-score per pathway (column-wise) for visualization
gsva_z <- scale(gsva_mat)
cat("  Heatmap matrix:", nrow(gsva_z), "rows x", ncol(gsva_z), "pathways\n")
cat("  Z-score range:", round(min(gsva_z, na.rm = TRUE), 2), "to",
    round(max(gsva_z, na.rm = TRUE), 2), "\n")

# Save scores
gsva_df <- as.data.frame(gsva_z) %>% tibble::rownames_to_column("celltype_age")
fwrite(gsva_df, file.path(output_dir, "multicelltype_pathway_scores.csv"))
cat("  Saved multicelltype_pathway_scores.csv\n")

# -----------------------------------------------------------------------------
# Step 5: Generate Figure 7 B/C heatmaps
# -----------------------------------------------------------------------------
cat("\nStep 5: Generating Figure 7 B/C dual heatmaps...\n")

gsva_t <- gsva_z  # rows = celltype_age, cols = pathways

# Row annotations
get_category <- function(ct_age) {
  ct <- gsub("_(Elderly|Young)$", "", ct_age)
  dplyr::case_when(
    grepl("Tcells|NK|Bcells|Plasma|Cycling", ct) & !grepl("Myeloid", ct) ~ "Lymphocyte",
    grepl("Macro|Mono|DC|Myeloid", ct) ~ "Myeloid",
    grepl("Epithelial|Cancer", ct) ~ "Epithelial",
    grepl("CAF|PVL|Endo", ct) ~ "Stromal",
    TRUE ~ "Other"
  )
}

row_categories <- data.frame(
  Category = sapply(rownames(gsva_t), get_category),
  row.names = rownames(gsva_t)
)

category_order <- c("Lymphocyte", "Myeloid", "Epithelial", "Stromal")
row_order <- rownames(gsva_t)[order(
  match(row_categories$Category, category_order),
  gsub("_(Elderly|Young)$", "", rownames(gsva_t)),
  grepl("_Elderly$", rownames(gsva_t))
)]
gsva_t <- gsva_t[row_order, ]
row_categories <- row_categories[row_order, , drop = FALSE]

# Rename for display
rownames(gsva_t) <- gsub("_Young$", "_Younger", rownames(gsva_t))
rownames(gsva_t) <- gsub("_Elderly$", "_Older", rownames(gsva_t))
rownames(row_categories) <- rownames(gsva_t)

# Select pathways
hallmark_cols <- grep("^(HALLMARK_|LI_ESTROGENE_)", colnames(gsva_t), value = TRUE)
biocarta_cols <- grep("^BIOCARTA_", colnames(gsva_t), value = TRUE)

if (pathway_mode == "divergent") {
  cat("  Selecting divergent pathways...\n")
  select_divergent <- function(mat, pw_cols, n_sel) {
    younger_rows <- grep("_Younger$", rownames(mat))
    older_rows <- grep("_Older$", rownames(mat))
    divergence <- sapply(pw_cols, function(pw) {
      abs(mean(mat[older_rows, pw], na.rm = TRUE) - mean(mat[younger_rows, pw], na.rm = TRUE))
    })
    sorted <- names(sort(divergence, decreasing = TRUE))
    selected <- head(sorted, n_sel)
    cat(sprintf("    Top 5 of %d:\n", length(selected)))
    for (i in 1:min(5, length(selected))) {
      cat(sprintf("      %d. %s (div=%.2f)\n", i, selected[i], divergence[selected[i]]))
    }
    selected
  }
  hallmark_keep <- select_divergent(gsva_t, hallmark_cols, n_pathways)
  biocarta_keep <- select_divergent(gsva_t, biocarta_cols, n_pathways)
} else {
  hallmark_curated <- c(
    "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS",
    "LI_ESTROGENE_EARLY_E2_RESPONSE_UP", "LI_ESTROGENE_LATE_E2_RESPONSE_UP"
  )
  biocarta_curated <- c(
    "BIOCARTA_INFLAM_PATHWAY", "BIOCARTA_IL6_PATHWAY", "BIOCARTA_IL2_PATHWAY",
    "BIOCARTA_NFKB_PATHWAY", "BIOCARTA_TNFR1_PATHWAY", "BIOCARTA_DEATH_PATHWAY",
    "BIOCARTA_FAS_PATHWAY", "BIOCARTA_CASPASE_PATHWAY", "BIOCARTA_P53_PATHWAY",
    "BIOCARTA_CELLCYCLE_PATHWAY", "BIOCARTA_G1_PATHWAY", "BIOCARTA_G2_PATHWAY"
  )
  hallmark_keep <- intersect(hallmark_curated, hallmark_cols)
  biocarta_keep <- intersect(biocarta_curated, biocarta_cols)

  hallmark_missing <- setdiff(hallmark_curated, hallmark_cols)
  biocarta_missing <- setdiff(biocarta_curated, biocarta_cols)
  if (length(hallmark_missing) > 0) warning("Missing HALLMARK: ", paste(hallmark_missing, collapse = ", "))
  if (length(biocarta_missing) > 0) warning("Missing BIOCARTA: ", paste(biocarta_missing, collapse = ", "))
  cat("  Using", length(hallmark_keep), "HALLMARK and", length(biocarta_keep), "BIOCARTA curated pathways\n")
}

stopifnot(length(hallmark_keep) > 0, length(biocarta_keep) > 0)

gsva_hallmark <- gsva_t[, hallmark_keep, drop = FALSE]
gsva_biocarta <- gsva_t[, biocarta_keep, drop = FALSE]

# Add significance stars from stats_df
# Build a lookup: celltype_pathway -> padj
sig_lookup <- setNames(stats_df$padj, paste0(stats_df$celltype, "::", stats_df$pathway))

# Function to create significance annotation matrix
make_sig_matrix <- function(mat, pw_names_original) {
  sig_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
  for (i in 1:nrow(mat)) {
    rn <- rownames(mat)[i]
    ct <- gsub("_(Younger|Older)$", "", rn)
    for (j in 1:ncol(mat)) {
      key <- paste0(ct, "::", pw_names_original[j])
      padj <- sig_lookup[key]
      if (!is.na(padj) && padj < 0.05) sig_mat[i, j] <- "*"
      if (!is.na(padj) && padj < 0.01) sig_mat[i, j] <- "**"
    }
  }
  sig_mat
}

hallmark_sig <- make_sig_matrix(gsva_hallmark, hallmark_keep)
biocarta_sig <- make_sig_matrix(gsva_biocarta, biocarta_keep)

# Clean column names for display
colnames(gsva_hallmark) <- gsub("^(HALLMARK_|LI_ESTROGENE_)", "", colnames(gsva_hallmark))
colnames(gsva_hallmark) <- gsub("_", " ", colnames(gsva_hallmark))
colnames(gsva_biocarta) <- gsub("^BIOCARTA_", "", colnames(gsva_biocarta))
colnames(gsva_biocarta) <- gsub("_", " ", colnames(gsva_biocarta))

ann_colors <- list(
  Category = c(Lymphocyte = "#4DAF4A", Myeloid = "#E41A1C", Epithelial = "#377EB8", Stromal = "#984EA3")
)

hallmark_col_fun <- colorRamp2(c(-3, -1.5, 0, 1.5, 3), c("#3E5CA8", "#8BC5E4", "#EDF3F1", "#F6BD6A", "#E75321"))
biocarta_col_fun <- colorRamp2(c(-4, -2, 0, 2, 4), c("#3E5CA8", "#8BC5E4", "#EDF3F1", "#F6BD6A", "#E75321"))

row_split <- factor(row_categories$Category, levels = category_order)

left_anno <- rowAnnotation(
  Category = row_categories$Category,
  col = ann_colors,
  show_annotation_name = FALSE,
  width = unit(3, "mm"),
  border = FALSE
)

make_heatmap <- function(mat, sig_mat, title_text, col_fun, range_vals) {
  Heatmap(
    mat, name = "z-score", col = col_fun,
    cluster_rows = FALSE, cluster_columns = TRUE,
    show_row_dend = FALSE,
    row_split = row_split, row_gap = unit(2, "mm"),
    row_title_side = "left", row_title_rot = 0,
    row_title_gp = gpar(fontsize = 12, fontface = "plain", fontfamily = "Arial"),
    left_annotation = left_anno,
    show_row_names = TRUE, row_names_side = "right",
    row_names_gp = gpar(fontsize = 11, fontfamily = "Arial"),
    show_column_names = TRUE, column_names_rot = 45,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 12, fontfamily = "Arial"),
    column_title = title_text,
    column_title_gp = gpar(fontsize = 18, fontface = "plain", fontfamily = "Arial"),
    column_dend_height = unit(10, "mm"),
    border = FALSE,
    rect_gp = gpar(col = "#3B3B3B", lwd = 0.5),
    # Overlay significance stars
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (sig_mat[i, j] != "") {
        grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 10, col = "black", fontfamily = "Arial"))
      }
    },
    width = unit(55, "mm"),
    height = unit(110, "mm"),
    heatmap_legend_param = list(
      direction = "horizontal",
      title = "Pathway activity\n(z-score)",
      at = c(range_vals[1], 0, range_vals[2]),
      legend_width = unit(30, "mm"),
      title_gp = gpar(fontsize = 12, fontfamily = "Arial"),
      labels_gp = gpar(fontsize = 11, fontfamily = "Arial")
    )
  )
}

ht_hallmark <- make_heatmap(gsva_hallmark, hallmark_sig, "Hallmark", hallmark_col_fun, c(-3, 3))
ht_biocarta <- make_heatmap(gsva_biocarta, biocarta_sig, "BIOCARTA", biocarta_col_fun, c(-4, 4))

mode_suffix <- ifelse(pathway_mode == "divergent", "_divergent", "")

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
