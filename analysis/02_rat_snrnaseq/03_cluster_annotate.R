#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/03_cluster_annotate.R
# UMAP, clustering, and cell type annotation for rat snRNA-seq
#
# Inputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_integrated.rds
#
# Outputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_annotated.rds
#   - analysis/02_rat_snrnaseq/outputs/clustering_plots.pdf
#   - analysis/02_rat_snrnaseq/outputs/marker_heatmap.pdf

# Set random seed for reproducibility
set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
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
fig_dir <- file.path(project_root, "figures/by_analysis/rat_snrnaseq")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Rat snRNA-seq Clustering and Annotation ===\n")
cat("Random seed: 12345\n")
cat("Project root:", project_root, "\n")
cat("Output directory:", output_dir, "\n\n")

# Input file
input_rds <- file.path(output_dir, "seurat_integrated.rds")

# Verify input file exists
cat("Checking input files...\n")
check_file_exists(input_rds, "Integrated Seurat object")
cat("\n")

# -----------------------------------------------------------------------------
# Step 1: Load Integrated Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading integrated Seurat object...\n")
seurat_obj <- readRDS(input_rds)
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Features:", nrow(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Run UMAP (if not already present or to ensure reproducibility)
# -----------------------------------------------------------------------------
cat("\nStep 2: Running UMAP...\n")

# Use Harmony embeddings for UMAP
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  seed.use = 12345,
  verbose = FALSE
)

cat("  UMAP complete\n")

# -----------------------------------------------------------------------------
# Step 3: Find Neighbors and Clusters
# -----------------------------------------------------------------------------
cat("\nStep 3: Finding neighbors and clustering...\n")

seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

# Run clustering at multiple resolutions
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (res in resolutions) {
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = res,
    random.seed = 12345,
    verbose = FALSE
  )
}

# Set default resolution
default_res <- 0.4
Idents(seurat_obj) <- paste0("SCT_snn_res.", default_res)
seurat_obj$seurat_clusters <- seurat_obj@meta.data[[paste0("SCT_snn_res.", default_res)]]

cat("  Clustering complete at resolutions:", paste(resolutions, collapse = ", "), "\n")
cat("  Default resolution:", default_res, "\n")
cat("  Number of clusters:", length(unique(seurat_obj$seurat_clusters)), "\n")

# -----------------------------------------------------------------------------
# Step 4: Define Cell Type Markers
# -----------------------------------------------------------------------------
cat("\nStep 4: Defining cell type markers...\n")

# Canonical markers for rat mammary tissue
# Note: Adjust gene names for rat nomenclature if needed
cell_type_markers <- list(
  # Epithelial cells
  Luminal = c("Krt8", "Krt18", "Krt19", "Epcam", "Cd24a"),
  Basal = c("Krt14", "Krt5", "Krt17", "Trp63", "Acta2"),

  # Immune cells
  T_cells = c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a"),
  B_cells = c("Cd19", "Cd79a", "Cd79b", "Ms4a1"),
  Macrophages = c("Cd68", "Adgre1", "Csf1r", "Mrc1"),

  # Stromal cells
  Fibroblasts = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Pdgfra"),
  Endothelial = c("Pecam1", "Cdh5", "Vwf", "Kdr"),

  # Adipocytes
  Adipocytes = c("Adipoq", "Pparg", "Fabp4", "Lep")
)

cat("  Defined markers for", length(cell_type_markers), "cell types\n")

# -----------------------------------------------------------------------------
# Step 5: Calculate Marker Expression per Cluster
# -----------------------------------------------------------------------------
cat("\nStep 5: Calculating marker expression per cluster...\n")

# Switch to RNA assay for marker expression
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

# Get average expression of markers per cluster
all_markers <- unlist(cell_type_markers)
available_markers <- all_markers[all_markers %in% rownames(seurat_obj)]
cat("  Available markers:", length(available_markers), "/", length(all_markers), "\n")

avg_exp <- AverageExpression(
  seurat_obj,
  features = available_markers,
  group.by = "seurat_clusters"
)$RNA

# -----------------------------------------------------------------------------
# Step 6: Assign Cell Types Based on Marker Expression
# -----------------------------------------------------------------------------
cat("\nStep 6: Assigning cell types based on marker expression...\n")

# Calculate scores for each cell type
cluster_ids <- colnames(avg_exp)
cell_type_assignments <- data.frame(cluster = cluster_ids, stringsAsFactors = FALSE)

for (ct in names(cell_type_markers)) {
  markers <- cell_type_markers[[ct]]
  markers_present <- markers[markers %in% rownames(avg_exp)]
  if (length(markers_present) > 0) {
    cell_type_assignments[[ct]] <- colMeans(avg_exp[markers_present, , drop = FALSE])
  } else {
    cell_type_assignments[[ct]] <- 0
  }
}

# Assign cell type based on highest scoring marker set
score_cols <- names(cell_type_markers)
cell_type_assignments$assigned_type <- apply(
  cell_type_assignments[, score_cols],
  1,
  function(x) score_cols[which.max(x)]
)

cat("  Cell type assignments:\n")
print(table(cell_type_assignments$assigned_type))

# Map cluster to cell type in Seurat object
cluster_to_celltype <- setNames(
  cell_type_assignments$assigned_type,
  cell_type_assignments$cluster
)

seurat_obj$CellTypeByMarker_RatsnRNAseq <- as.character(
  cluster_to_celltype[as.character(seurat_obj$seurat_clusters)]
)

# Also create broader macrophage/T cell subtype column for compatibility
seurat_obj$CellTypeMacroTcell_RatsnRNAseq <- seurat_obj$CellTypeByMarker_RatsnRNAseq

# -----------------------------------------------------------------------------
# Step 7: Generate Clustering and Annotation Plots
# -----------------------------------------------------------------------------
cat("\nStep 7: Generating clustering and annotation plots...\n")

pdf(file.path(output_dir, "clustering_plots.pdf"), width = 14, height = 10)

# UMAP colored by cluster
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle(paste0("UMAP by Cluster (resolution=", default_res, ")")) +
  theme_bw()
print(p1)

# UMAP colored by cell type
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "CellTypeByMarker_RatsnRNAseq", label = TRUE) +
  ggtitle("UMAP by Cell Type") +
  theme_bw()
print(p2)

# UMAP split by age group
p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "CellTypeByMarker_RatsnRNAseq",
              split.by = "AgeGroup") +
  ggtitle("UMAP by Cell Type - Split by Age Group") +
  theme_bw()
print(p3)

# UMAP at different resolutions
for (res in resolutions) {
  col_name <- paste0("SCT_snn_res.", res)
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = col_name, label = TRUE) +
    ggtitle(paste0("UMAP at resolution=", res)) +
    theme_bw()
  print(p)
}

dev.off()

cat("  Clustering plots saved to:", file.path(output_dir, "clustering_plots.pdf"), "\n")

# Marker heatmap
pdf(file.path(output_dir, "marker_heatmap.pdf"), width = 12, height = 10)

# Scale the average expression for heatmap
avg_exp_scaled <- t(scale(t(avg_exp)))

# Create annotation for clusters
cluster_annotation <- data.frame(
  CellType = cluster_to_celltype[colnames(avg_exp_scaled)],
  row.names = colnames(avg_exp_scaled)
)

pheatmap(
  avg_exp_scaled,
  main = "Cell Type Marker Expression by Cluster",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = cluster_annotation,
  fontsize = 8
)

dev.off()

cat("  Marker heatmap saved to:", file.path(output_dir, "marker_heatmap.pdf"), "\n")

# -----------------------------------------------------------------------------
# Step 8: Find Cluster Markers (for verification)
# -----------------------------------------------------------------------------
cat("\nStep 8: Finding cluster markers...\n")

Idents(seurat_obj) <- "seurat_clusters"
cluster_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  test.use = "wilcox",
  verbose = FALSE
)

# Save top markers per cluster
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

write.csv(top_markers, file.path(output_dir, "cluster_markers_top10.csv"), row.names = FALSE)
cat("  Top 10 markers per cluster saved\n")

# -----------------------------------------------------------------------------
# Step 9: Save Annotated Object
# -----------------------------------------------------------------------------
cat("\nStep 9: Saving annotated object...\n")
saveRDS(seurat_obj, file.path(output_dir, "seurat_annotated.rds"))

# Also save cell type assignments
write.csv(cell_type_assignments, file.path(output_dir, "cell_type_assignments.csv"), row.names = FALSE)

# Summary table
summary_table <- seurat_obj@meta.data %>%
  group_by(CellTypeByMarker_RatsnRNAseq, AgeGroup) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = AgeGroup, values_from = n_cells, values_fill = 0)

write.csv(summary_table, file.path(output_dir, "cell_counts_by_type_age.csv"), row.names = FALSE)
cat("  Cell count summary saved\n")

# Memory cleanup
gc()

cat("\n=== Clustering and annotation complete ===\n")
cat("Output:", file.path(output_dir, "seurat_annotated.rds"), "\n")
cat("\nNOTE: Cell type assignments are automated based on marker expression.\n")
cat("Manual review and refinement may be needed based on biological knowledge.\n")
