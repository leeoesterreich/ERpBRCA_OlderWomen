#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/03_cell_type_annotation.R
# Map Xu et al. 2024 cell type annotations to pipeline-standard names
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_young_midage_elderly.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#   - analysis/04_human_scrnaseq/outputs/umap_celltypes.pdf
#   - analysis/04_human_scrnaseq/outputs/metadata_annotated.txt
#   - analysis/04_human_scrnaseq/figures/umap_celltypes.png

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(patchwork)
})

# Set Arial as default font for all plots
library(showtext)
font_add("Arial", "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Arial.ttf")
showtext_auto()
theme_set(theme_bw(base_size = 14, base_family = "Arial"))

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

cat("=== Cell Type Annotation (Xu et al. 2024 → Pipeline Names) ===\n")

# -----------------------------------------------------------------------------
# Cell type mapping: Xu annotation → pipeline standard names
# -----------------------------------------------------------------------------
CELLTYPE_MAP <- c(
  "Cancer Epithelial Cells"       = "CancerEpithelial",
  "CD4+ T Cells"                  = "TcellsCD4",
  "CD8+ T Cells"                  = "TcellsCD8",
  "Regulatory T Cells"            = "Tregs",
  "NK Cells"                      = "NKcells",
  "B Cells"                       = "Bcells",
  "Plasma Cells"                  = "Plasmablasts",
  "Macrophages"                   = "Macrophage",
  "Monocytes"                     = "Monocyte",
  "Dendritic Cells"               = "DCs",
  "Fibroblasts"                   = "CAFs",
  "Endothelial Cells"             = "Endothelial",
  "Perivascular-like (PVL) Cells" = "PVL",
  "Epithelial Cells"              = "NormalEpithelial",
  "Myoepithelial Cells"           = "Myoepithelial",
  "MDSCs"                         = "MDSCs",
  "Mast Cells"                    = "MastCells",
  "Neutrophils"                   = "Neutrophils"
)

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_young_midage_elderly.rds"))
cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Map cell type annotations
# -----------------------------------------------------------------------------
cat("\nStep 2: Mapping cell types...\n")

# Get original Xu annotations
xu_types <- seurat_obj$CellTypeMajor
cat("  Original Xu cell types:\n")
print(table(xu_types, useNA = "ifany"))

# Map to pipeline names
mapped_types <- CELLTYPE_MAP[xu_types]

# Check for unmapped types
unmapped <- xu_types[is.na(mapped_types)]
if (length(unmapped) > 0) {
  cat("\n  WARNING: Unmapped cell types found:\n")
  print(table(unmapped))
  cat("  These cells will be labeled 'Unknown'\n")
  mapped_types[is.na(mapped_types)] <- "Unknown"
}

seurat_obj$CellTypeAnnot <- as.character(mapped_types)
seurat_obj$CellTypeAnnotSH <- seurat_obj$CellTypeAnnot

cat("\n  Mapped cell types:\n")
print(table(seurat_obj$CellTypeAnnot))

# -----------------------------------------------------------------------------
# Step 3: Set identity and visualize
# -----------------------------------------------------------------------------
cat("\nStep 3: Visualizing...\n")

Idents(seurat_obj) <- seurat_obj$CellTypeAnnot

p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4) +
  ggtitle("Cell Types") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "AgeGroup",
              cols = c(Elderly = "#D55E00", MidAge = "#009E73", Young = "#56B4E9")) +
  ggtitle("Age Groups") +
  theme_classic(base_size = 14)

combined <- wrap_plots(list(p1, p2), ncol = 2)

pdf(file.path(output_dir, "umap_celltypes.pdf"), width = 14, height = 6)
print(combined)
dev.off()

ggsave(file.path(figures_dir, "umap_celltypes.png"), combined,
       width = 14, height = 6, dpi = 300, bg = "white")
ggsave(file.path(figures_dir, "umap_celltypes.svg"), combined,
       width = 14, height = 6)  # SVG for vector assembly

# Feature plots for key markers
markers <- c("EPCAM", "KRT19", "CD68", "CD3D", "MS4A1", "PECAM1",
             "CTLA4", "PDCD1", "CCL2", "TGFB1", "MRC1", "CD163")
markers_present <- markers[markers %in% rownames(seurat_obj)]

if (length(markers_present) > 0) {
  pdf(file.path(output_dir, "feature_markers.pdf"), width = 16, height = 12)
  print(FeaturePlot(seurat_obj, features = markers_present, ncol = 4))
  dev.off()
}

# -----------------------------------------------------------------------------
# Step 4: Save metadata
# -----------------------------------------------------------------------------
cat("\nStep 4: Saving...\n")

saveRDS(seurat_obj, file.path(output_dir, "seurat_annotated.rds"))

metadata <- seurat_obj@meta.data %>%
  tibble::rownames_to_column("CellID") %>%
  select(CellID, orig.ident, AgeGroup, CellTypeAnnot, CellTypeAnnotSH, Dataset)

fwrite(metadata, file.path(output_dir, "metadata_annotated.txt"),
       sep = "\t", quote = FALSE)

cat("  Saved seurat_annotated.rds:", ncol(seurat_obj), "cells\n")
cat("  Saved metadata_annotated.txt\n")

cat("\n=== Annotation complete ===\n")
