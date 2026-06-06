#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/10_multicelltype_pathway.R
# Figure 2E - HALLMARK pathway expression across immune cell subtypes
# in tumors from younger vs older rats (snRNA-seq)
#
# Mirrors analysis/04_human_scrnaseq/15_multicelltype_pathway.R but for rat.
#
# Inputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_annotated_sctype.rds (Seurat obj with
#     `sctype_celltype` metadata column + `AgeGroup` = Young|Aged)
#
# Outputs:
#   - analysis/02_rat_snrnaseq/outputs/rat_multicelltype_pathway_scores.csv
#   - analysis/02_rat_snrnaseq/outputs/rat_multicelltype_pathway_stats.csv
#   - analysis/02_rat_snrnaseq/figures/fig2e_rat_pathway_heatmap.{svg,pdf,png}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(msigdbr)
  library(GSVA)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(showtext)
})

tryCatch(
  font_add("Arial", "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"),
  error = function(e) message("Arial not found; using default sans font")
)
showtext_auto()

cat("=== Rat snRNA-seq multicelltype pathway analysis (Fig 2E) ===\n")

output_dir <- "analysis/02_rat_snrnaseq/outputs"
figures_dir <- "analysis/02_rat_snrnaseq/figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Step 1: Load annotated Seurat object
# -----------------------------------------------------------------------------
cat("\nStep 1: Loading rat annotated Seurat object...\n")
seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated_sctype.rds"))
cat("  Cells:", ncol(seurat_obj), "  Genes:", nrow(seurat_obj), "\n")

# Pull metadata
md <- seurat_obj@meta.data
stopifnot(all(c("sctype_celltype", "AgeGroup", "orig.ident") %in% colnames(md)))

# Restrict to immune + key tumor cell types relevant to Fig 2E
keep_types <- c(
  "CancerEpithelial",
  "CD4Tcell", "CD8Tcell", "NaiveTcell", "Treg",
  "NKcell",
  "Macrophage", "M1_Macrophage", "M2_Macrophage", "Monocyte", "DendriticCell",
  "Endothelial", "Fibroblast"
)
md$keep <- md$sctype_celltype %in% keep_types
seurat_obj <- subset(seurat_obj, cells = rownames(md)[md$keep])
cat("  Cells after filter:", ncol(seurat_obj), "\n")

# Normalize Age label to match human convention (Young / Elderly)
seurat_obj$AgeGroup <- ifelse(seurat_obj$AgeGroup == "Aged", "Elderly", "Young")

# -----------------------------------------------------------------------------
# Step 2: Load HALLMARK gene sets (rat)
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading rat HALLMARK gene sets...\n")
hallmark_df <- tryCatch(
  msigdbr(species = "Rattus norvegicus", collection = "H"),
  error = function(e) msigdbr(species = "Rattus norvegicus", category = "H")
)
gs_col <- if ("gs_name" %in% colnames(hallmark_df)) "gs_name" else "gs_collection"
sym_col <- if ("gene_symbol" %in% colnames(hallmark_df)) "gene_symbol" else "gs_gene_symbol"
hallmark_list <- split(hallmark_df[[sym_col]], hallmark_df[[gs_col]])
cat("  Loaded", length(hallmark_list), "HALLMARK gene sets\n")

# -----------------------------------------------------------------------------
# Step 3: Pseudo-bulk per sample x cell type (orig.ident is per-rat sample)
# -----------------------------------------------------------------------------
cat("\nStep 3: Building pseudo-bulk matrix (sample x celltype x age)...\n")
DefaultAssay(seurat_obj) <- "RNA"
# Seurat 5: counts may be split into per-sample layers; join into one
if (inherits(seurat_obj[["RNA"]], "Assay5")) {
  cat("  Joining Seurat 5 multi-sample counts layers...\n")
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
}
counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

md <- seurat_obj@meta.data
md$ct_age <- paste0(md$sctype_celltype, "_", md$AgeGroup)
md$ct_age_sample <- paste0(md$ct_age, "::", md$orig.ident)

# Aggregate to pseudo-bulk: sum counts per (cell_type, age, sample)
groups <- unique(md$ct_age_sample)
cat("  Pseudo-bulk groups:", length(groups), "\n")

# Use rowSums per group; skip groups with too few cells
min_cells <- 5
pb_list <- list()
group_meta <- list()
for (g in groups) {
  cells <- rownames(md)[md$ct_age_sample == g]
  if (length(cells) < min_cells) next
  pb_list[[g]] <- rowSums(counts[, cells, drop = FALSE])
  group_meta[[g]] <- list(
    celltype = md$sctype_celltype[match(cells[1], rownames(md))],
    age      = md$AgeGroup[match(cells[1], rownames(md))],
    sample   = md$orig.ident[match(cells[1], rownames(md))],
    n_cells  = length(cells)
  )
}
pb_mat <- do.call(cbind, pb_list)
colnames(pb_mat) <- names(pb_list)
cat("  Pseudo-bulk matrix:", nrow(pb_mat), "genes x", ncol(pb_mat), "groups\n")

# CPM + log2
lib_size <- colSums(pb_mat)
cpm <- sweep(pb_mat, 2, lib_size, FUN = "/") * 1e6
log_cpm <- log2(cpm + 1)

# -----------------------------------------------------------------------------
# Step 4: GSVA per cell type (samples as columns within each cell type group)
# We instead run GSVA across ALL groups in one shot, then split for stats
# -----------------------------------------------------------------------------
cat("\nStep 4: Running GSVA...\n")
gsva_par <- tryCatch(
  gsvaParam(exprData = as.matrix(log_cpm), geneSets = hallmark_list, kcdf = "Gaussian"),
  error = function(e) NULL
)
gsva_mat <- if (!is.null(gsva_par)) GSVA::gsva(gsva_par) else GSVA::gsva(as.matrix(log_cpm), hallmark_list, kcdf = "Gaussian")
cat("  GSVA matrix:", nrow(gsva_mat), "pathways x", ncol(gsva_mat), "groups\n")

# -----------------------------------------------------------------------------
# Step 5: Aggregate to celltype_age (mean across samples within group)
# Then z-score across (celltype_age) for each pathway
# -----------------------------------------------------------------------------
group_df <- do.call(rbind, lapply(names(group_meta), function(g) {
  data.frame(group = g,
             celltype = group_meta[[g]]$celltype,
             age = group_meta[[g]]$age,
             sample = group_meta[[g]]$sample,
             n_cells = group_meta[[g]]$n_cells,
             stringsAsFactors = FALSE)
}))

# Mean per (celltype, age)
gsva_t <- t(gsva_mat)  # rows = groups, cols = pathways
gsva_df <- as.data.frame(gsva_t)
gsva_df$group <- rownames(gsva_df)
gsva_long <- gsva_df %>%
  tidyr::pivot_longer(-group, names_to = "pathway", values_to = "score") %>%
  dplyr::left_join(group_df, by = "group") %>%
  dplyr::group_by(celltype, age, pathway) %>%
  dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(celltype_age = paste0(celltype, "_", age))

# Pivot to celltype_age x pathway
gsva_ca <- gsva_long %>%
  dplyr::select(celltype_age, pathway, mean_score) %>%
  tidyr::pivot_wider(names_from = pathway, values_from = mean_score) %>%
  as.data.frame()
rownames(gsva_ca) <- gsva_ca$celltype_age
gsva_ca$celltype_age <- NULL
gsva_ca_mat <- as.matrix(gsva_ca)

# z-score per pathway (column) across cell types/ages
gsva_z <- scale(gsva_ca_mat, center = TRUE, scale = TRUE)

# Save raw scores + z-scored
fwrite(
  data.frame(celltype_age = rownames(gsva_ca_mat), gsva_ca_mat),
  file.path(output_dir, "rat_multicelltype_pathway_scores_raw.csv")
)
fwrite(
  data.frame(celltype_age = rownames(gsva_z), gsva_z),
  file.path(output_dir, "rat_multicelltype_pathway_scores.csv")
)
cat("  Saved rat_multicelltype_pathway_scores.csv\n")

# -----------------------------------------------------------------------------
# Step 6: Per-pathway, per-cell-type Welch t-test Young vs Elderly across samples
# -----------------------------------------------------------------------------
cat("\nStep 5: Per-celltype Welch t-test (Young vs Elderly)...\n")
stats_list <- list()
for (ct in unique(group_df$celltype)) {
  ct_groups <- group_df$group[group_df$celltype == ct]
  if (length(ct_groups) < 2) next
  ct_age <- group_df$age[group_df$celltype == ct]
  if (sum(ct_age == "Young") < 2 || sum(ct_age == "Elderly") < 2) next
  for (pw in rownames(gsva_mat)) {
    young_vals <- gsva_mat[pw, ct_groups[ct_age == "Young"]]
    elderly_vals <- gsva_mat[pw, ct_groups[ct_age == "Elderly"]]
    if (length(young_vals) < 2 || length(elderly_vals) < 2) next
    p <- tryCatch(t.test(elderly_vals, young_vals)$p.value, error = function(e) NA)
    stats_list[[length(stats_list) + 1]] <- data.frame(
      celltype = ct, pathway = pw,
      young_mean = mean(young_vals), elderly_mean = mean(elderly_vals),
      diff = mean(elderly_vals) - mean(young_vals),
      n_young = length(young_vals), n_elderly = length(elderly_vals),
      pvalue = p
    )
  }
}
stats_df <- do.call(rbind, stats_list)
stats_df$padj <- p.adjust(stats_df$pvalue, method = "BH")
stats_df$significant <- !is.na(stats_df$padj) & stats_df$padj < 0.05
fwrite(stats_df, file.path(output_dir, "rat_multicelltype_pathway_stats.csv"))
cat("  Saved rat_multicelltype_pathway_stats.csv\n")

# -----------------------------------------------------------------------------
# Step 7: Render Fig 2E heatmap (rat HALLMARK x CellType_Age)
# -----------------------------------------------------------------------------
cat("\nStep 6: Rendering Fig 2E heatmap...\n")
mat <- t(gsva_z)  # rows = pathways, cols = celltype_age

get_category <- function(ct) {
  dplyr::case_when(
    grepl("CD4Tcell|CD8Tcell|NaiveTcell|NKcell|Treg", ct) ~ "Lymphocyte",
    grepl("Macrophage|Monocyte|DendriticCell", ct) ~ "Myeloid",
    grepl("CancerEpithelial|Myoepithelial", ct) ~ "Epithelial",
    grepl("Endothelial|Fibroblast", ct) ~ "Stromal",
    TRUE ~ "Other"
  )
}

col_celltypes <- sub("_(Young|Elderly)$", "", colnames(mat))
col_categories <- sapply(col_celltypes, get_category)
category_order <- c("Lymphocyte", "Myeloid", "Epithelial", "Stromal", "Other")

col_order <- order(match(col_categories, category_order),
                   col_celltypes,
                   grepl("_Elderly$", colnames(mat)))
mat <- mat[, col_order, drop = FALSE]
col_categories <- col_categories[col_order]

# Display: Young -> Younger, Elderly -> Older
display_cols <- gsub("_Young$", "_Younger", colnames(mat))
display_cols <- gsub("_Elderly$", "_Older", display_cols)
colnames(mat) <- display_cols

# Strip HALLMARK_ prefix, _ -> space
display_rows <- gsub("^HALLMARK_", "", rownames(mat))
display_rows <- gsub("_", " ", display_rows)
rownames(mat) <- display_rows

# Significance stars
sig_lookup <- setNames(stats_df$padj, paste0(stats_df$celltype, "::", stats_df$pathway))
sig_mat <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
for (j in seq_len(ncol(mat))) {
  ct_age <- colnames(mat)[j]
  ct <- gsub("_(Younger|Older)$", "", ct_age)
  for (i in seq_len(nrow(mat))) {
    pw <- paste0("HALLMARK_", gsub(" ", "_", rownames(mat)[i]))
    key <- paste0(ct, "::", pw)
    padj <- sig_lookup[key]
    if (!is.na(padj) && padj < 0.05) sig_mat[i, j] <- "*"
    if (!is.na(padj) && padj < 0.01) sig_mat[i, j] <- "**"
  }
}

ann_colors <- list(
  Category = c(Lymphocyte = "#009E73", Myeloid = "#D55E00",
               Epithelial = "#0072B2", Stromal = "#CC79A7", Other = "#000000")
)
col_split <- factor(col_categories, levels = category_order)
hallmark_col_fun <- colorRamp2(
  c(-3, -1.5, 0, 1.5, 3),
  c("#3E5CA8", "#8BC5E4", "#EDF3F1", "#F6BD6A", "#E75321")
)

bottom_anno <- HeatmapAnnotation(
  Category = col_categories,
  col = ann_colors,
  show_annotation_name = FALSE,
  height = unit(3, "mm"),
  border = FALSE,
  which = "column"
)

ht <- Heatmap(
  mat,
  name = "Pathway activity",
  col = hallmark_col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_split, column_gap = unit(2, "mm"),
  column_title_gp = gpar(fontsize = 11, fontfamily = "Arial"),
  bottom_annotation = bottom_anno,
  show_row_names = TRUE, row_names_side = "right",
  row_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  show_column_names = TRUE, column_names_rot = 45,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 9, fontfamily = "Arial"),
  border = FALSE,
  rect_gp = gpar(col = "#3B3B3B", lwd = 0.4),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] != "") {
      grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 8, col = "black", fontfamily = "Arial"))
    }
  },
  width = unit(140, "mm"),
  height = unit(170, "mm"),
  heatmap_legend_param = list(
    direction = "horizontal",
    title = "Pathway activity\n(z-score)",
    at = c(-3, 0, 3),
    legend_width = unit(30, "mm"),
    title_gp = gpar(fontsize = 10, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 9, fontfamily = "Arial")
  )
)

svg_path <- file.path(figures_dir, "fig2e_rat_pathway_heatmap.svg")
pdf_path <- file.path(figures_dir, "fig2e_rat_pathway_heatmap.pdf")
png_path <- file.path(figures_dir, "fig2e_rat_pathway_heatmap.png")

svg(svg_path, width = 9.5, height = 9.0)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
invisible(dev.off())
cat("Saved:", svg_path, "\n")

pdf(pdf_path, width = 9.5, height = 9.0)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
invisible(dev.off())
cat("Saved:", pdf_path, "\n")

png(png_path, width = 9.5, height = 9.0, units = "in", res = 300)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
invisible(dev.off())
cat("Saved:", png_path, "\n")

cat("\n=== Rat 2E pipeline complete ===\n")
