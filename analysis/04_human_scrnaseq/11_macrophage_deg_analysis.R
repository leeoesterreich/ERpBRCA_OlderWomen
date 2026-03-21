#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/11_macrophage_deg_analysis.R
# Pseudobulk DEG analysis: macrophages from elderly vs young patients
# Aggregates cells per patient → DESeq2 (avoids pseudoreplication)
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/macrophage_seurat.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/macrophage_degs.csv
#   - analysis/04_human_scrnaseq/figures/fig7b_macrophage_volcano.png
#   - analysis/04_human_scrnaseq/figures/fig7b_macrophage_heatmap.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(DESeq2)
  library(pheatmap)
  library(ggrepel)
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
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Remove stale heatmap artifacts — prevents showing genes from prior runs
# when current run has 0 significant DEGs
stale_files <- list.files(figures_dir, pattern = "fig7b_macrophage_heatmap", full.names = TRUE)
if (length(stale_files) > 0) {
  cat("  Removing", length(stale_files), "stale heatmap files\n")
  file.remove(stale_files)
}

MIN_CELLS_PER_PSEUDOBULK <- 10

cat("=== Macrophage Pseudobulk DEG Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_file <- file.path(output_dir, "macrophage_seurat.rds")
if (!file.exists(seurat_file)) {
  stop("Run 10_load_macrophage_seurat.R first")
}

seurat_obj <- readRDS(seurat_file)
cat("  Cells:", ncol(seurat_obj), "\n")

cat("  Age group distribution:\n")
print(table(seurat_obj$AgeGroup))

# Subset to macrophages (in case object has both Macrophage + Monocyte)
macrophage_obj <- subset(seurat_obj, subset = CellTypeAnnotSH == "Macrophage")
cat("  Macrophages:", ncol(macrophage_obj), "cells\n")

# Filter to Elderly and Young
seurat_ey <- subset(macrophage_obj, subset = AgeGroup %in% c("Elderly", "Young"))
cat("  Macrophages Elderly vs Young:", ncol(seurat_ey), "cells\n")

cat("\n  Macrophages per age group:\n")
print(table(seurat_ey$AgeGroup))

cat("\n  Macrophages per patient:\n")
print(table(seurat_ey$orig.ident, seurat_ey$AgeGroup))

rm(seurat_obj, macrophage_obj); gc()

# -----------------------------------------------------------------------------
# Step 2: Pseudobulk aggregation
# -----------------------------------------------------------------------------
cat("\nStep 2: Aggregating to pseudobulk (per patient)...\n")

DefaultAssay(seurat_ey) <- "RNA"

# Remove patients with too few macrophages
cells_per_patient <- table(seurat_ey$orig.ident)
valid_patients <- names(cells_per_patient[cells_per_patient >= MIN_CELLS_PER_PSEUDOBULK])
cat("  Patients with >=", MIN_CELLS_PER_PSEUDOBULK, "macrophages:", length(valid_patients),
    "of", length(cells_per_patient), "\n")

seurat_ey <- subset(seurat_ey, cells = colnames(seurat_ey)[seurat_ey$orig.ident %in% valid_patients])

# Aggregate raw counts per patient
pseudo <- AggregateExpression(
  seurat_ey,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = "orig.ident"
)$RNA

# Seurat mangles column names: prepends "g" to IDs starting with digits,
# replaces "_" with "-". Reverse this to recover original Patient_IDs.
orig_ids <- unique(seurat_ey$orig.ident)
seurat_mangle <- function(x) {
  x <- gsub("_", "-", x)
  ifelse(grepl("^[0-9]", x), paste0("g", x), x)
}
pseudo_to_orig <- setNames(orig_ids, seurat_mangle(orig_ids))
colnames(pseudo) <- pseudo_to_orig[colnames(pseudo)]

# Warn about any unmatched columns
unmatched <- is.na(colnames(pseudo))
if (any(unmatched)) {
  cat("  WARNING: Could not map", sum(unmatched), "pseudobulk columns back to Patient_IDs\n")
  pseudo <- pseudo[, !unmatched, drop = FALSE]
}

# Build coldata
patient_age <- data.frame(
  orig.ident = unique(seurat_ey$orig.ident),
  AgeGroup = seurat_ey$AgeGroup[match(unique(seurat_ey$orig.ident), seurat_ey$orig.ident)],
  stringsAsFactors = FALSE
)
rownames(patient_age) <- patient_age$orig.ident
patient_age$orig.ident <- NULL
patient_age <- patient_age[patient_age$AgeGroup %in% c("Young", "Elderly"), , drop = FALSE]
patient_age$AgeGroup <- factor(patient_age$AgeGroup, levels = c("Young", "Elderly"))
pseudo <- pseudo[, rownames(patient_age)]

cat("  Pseudobulk matrix:", nrow(pseudo), "genes x", ncol(pseudo), "patients\n")
cat("  Young:", sum(patient_age$AgeGroup == "Young"),
    "patients, Elderly:", sum(patient_age$AgeGroup == "Elderly"), "patients\n")

# Filter low-count genes
keep <- rowSums(pseudo >= 10) >= 3
pseudo <- pseudo[keep, ]
cat("  Genes after filtering:", nrow(pseudo), "\n")

# -----------------------------------------------------------------------------
# Step 3: DESeq2
# -----------------------------------------------------------------------------
cat("\nStep 3: Running DESeq2...\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(pseudo),
  colData = patient_age,
  design = ~ AgeGroup
)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("AgeGroup", "Elderly", "Young"), alpha = 0.05)

degs <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

cat("  Total DEGs (FDR < 0.05):", sum(degs$padj < 0.05, na.rm = TRUE), "\n")
cat("  Upregulated in Elderly:", sum(degs$padj < 0.05 & degs$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Downregulated in Elderly:", sum(degs$padj < 0.05 & degs$log2FoldChange < 0, na.rm = TRUE), "\n")

fwrite(degs, file.path(output_dir, "macrophage_degs.csv"))
cat("\n  Saved macrophage_degs.csv\n")

cat("\n  Top 20 DEGs by FDR:\n")
print(head(degs, 20))

# -----------------------------------------------------------------------------
# Step 4: Volcano plot
# -----------------------------------------------------------------------------
cat("\nStep 4: Generating volcano plot...\n")

degs$significance <- case_when(
  degs$padj < 0.05 & degs$log2FoldChange > 0.5 ~ "Up",
  degs$padj < 0.05 & degs$log2FoldChange < -0.5 ~ "Down",
  TRUE ~ "NS"
)

top_genes <- degs %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(20) %>%
  pull(gene)

degs$label <- ifelse(degs$gene %in% top_genes, degs$gene, "")

volcano_plot <- ggplot(degs, aes(x = log2FoldChange, y = -log10(padj + 1e-300))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                  box.padding = 0.3, segment.color = "grey50") +
  scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey60")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(
    title = "Macrophage Pseudobulk DEGs: Elderly vs Young",
    subtitle = paste0("DESeq2 (n=", sum(patient_age$AgeGroup == "Young"), " Young, ",
                      sum(patient_age$AgeGroup == "Elderly"), " Elderly patients) | ",
                      sum(degs$padj < 0.05 & degs$log2FoldChange > 0, na.rm = TRUE), " up, ",
                      sum(degs$padj < 0.05 & degs$log2FoldChange < 0, na.rm = TRUE), " down"),
    x = "Log2 Fold Change (Elderly / Young)",
    y = "-Log10(FDR)",
    color = "Direction"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave(file.path(figures_dir, "fig7b_macrophage_volcano.png"),
       volcano_plot, width = 10, height = 8, dpi = 300)
cat("  Saved fig7b_macrophage_volcano.png\n")

# -----------------------------------------------------------------------------
# Step 5: Heatmap of top DEGs
# -----------------------------------------------------------------------------
cat("\nStep 5: Generating heatmap of top DEGs...\n")

top_degs <- degs %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(30) %>%
  pull(gene)

if (length(top_degs) > 0) {
  # Use normalized counts from DESeq2 for heatmap
  norm_counts <- counts(dds, normalized = TRUE)
  heatmap_mat <- norm_counts[top_degs[top_degs %in% rownames(norm_counts)], , drop = FALSE]

  ann_col <- data.frame(AgeGroup = patient_age$AgeGroup)
  rownames(ann_col) <- rownames(patient_age)
  ann_colors <- list(AgeGroup = c(Young = "#4DAF4A", Elderly = "#E41A1C"))

  png(file.path(figures_dir, "fig7b_macrophage_heatmap.png"),
      width = 10*300, height = 10*300, res = 300)
  pheatmap(
    as.matrix(heatmap_mat),
    scale = "row",
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    main = "Top Pseudobulk DEGs: Macrophages Elderly vs Young",
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    fontsize_row = 8
  )
  dev.off()
  cat("  Saved fig7b_macrophage_heatmap.png\n")
} else {
  cat("  No significant DEGs found for heatmap\n")
}

cat("\n=== Macrophage pseudobulk DEG analysis complete ===\n")
cat("Next: Run 12_macrophage_pathway_enrichment.R\n")
