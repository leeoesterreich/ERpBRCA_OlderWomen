#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/03_cell_type_annotation.R
# Cell type annotation using marker genes
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_young_midage_elderly.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#   - analysis/04_human_scrnaseq/outputs/umap_celltypes.pdf
#   - analysis/04_human_scrnaseq/outputs/metadata_annotated.txt

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(patchwork)
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

cat("=== Cell Type Annotation ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_young_midage_elderly.rds"))
cat("  Cells:", ncol(seurat_obj), "\n")

# Check if cell type annotations exist from original data
if ("CellTypeMajor" %in% colnames(seurat_obj@meta.data)) {
  cat("  Using existing cell type annotations\n")

  # Clean up cell type names (remove special characters)
  seurat_obj$CellTypeMajor <- gsub("-", "", seurat_obj$CellTypeMajor)

  # Create combined annotation
  # FIX: Ensure CellTypeMinor exists before using
  if ("CellTypeMinor" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$CellTypeAnnot <- ifelse(
      seurat_obj$CellTypeMajor %in% c("Myeloid", "Tcells"),
      seurat_obj$CellTypeMinor,
      seurat_obj$CellTypeMajor
    )
  } else {
    seurat_obj$CellTypeAnnot <- seurat_obj$CellTypeMajor
  }

  # Clean annotation names
  seurat_obj$CellTypeAnnot <- gsub("-|_|\\+", "", seurat_obj$CellTypeAnnot)

} else {
  cat("  No existing annotations, using cluster-based annotation\n")
  # Placeholder - would use marker-based annotation here
  seurat_obj$CellTypeAnnot <- paste0("Cluster_", seurat_obj$seurat_clusters)
}

cat("  Cell types:\n")
print(table(seurat_obj$CellTypeAnnot))

# -----------------------------------------------------------------------------
# Step 2: Set identity and visualize
# -----------------------------------------------------------------------------
cat("\nStep 2: Visualizing...\n")

Idents(seurat_obj) <- seurat_obj$CellTypeAnnot

p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4) +
  ggtitle("Cell Types") +
  theme(legend.position = "right")

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "AgeGroup") +
  ggtitle("Age Groups")

pdf(file.path(output_dir, "umap_celltypes.pdf"), width = 14, height = 6)
print(p1 + p2)
dev.off()

# Feature plots for key markers
markers <- c("EPCAM", "KRT19", "CD68", "CD3D", "MS4A1", "PECAM1")
markers_present <- markers[markers %in% rownames(seurat_obj)]

if (length(markers_present) > 0) {
  pdf(file.path(output_dir, "feature_markers.pdf"), width = 12, height = 8)
  print(FeaturePlot(seurat_obj, features = markers_present, ncol = 3))
  dev.off()
}

# -----------------------------------------------------------------------------
# Step 3: Save metadata
# -----------------------------------------------------------------------------
cat("\nStep 3: Saving...\n")

# Save annotated object
saveRDS(seurat_obj, file.path(output_dir, "seurat_annotated.rds"))

# Save metadata
metadata <- seurat_obj@meta.data %>%
  tibble::rownames_to_column("CellID") %>%
  select(CellID, orig.ident, AgeGroup, CellTypeAnnot)

fwrite(metadata, file.path(output_dir, "metadata_annotated.txt"),
       sep = "\t", quote = FALSE)

cat("\n=== Annotation complete ===\n")
