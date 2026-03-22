#!/usr/bin/env Rscript
# Spot check marker expression to validate scType annotations
# Compare with original Sanghoon marker expectations

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Set Arial as default font for all plots
library(showtext)
font_add("Arial", "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Arial.ttf")
showtext_auto()
theme_set(theme_bw(base_size = 14, base_family = "Arial"))

cat("=== Marker Expression Spot Check ===\n\n")

# Load scType cluster-level annotated object
obj_path <- "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/02_rat_snrnaseq/outputs/seurat_annotated_sctype_cluster.rds"
seurat_obj <- readRDS(obj_path)

cat("Loaded:", ncol(seurat_obj), "cells\n\n")

# Handle Seurat v5 layers
tryCatch({
  if ("JoinLayers" %in% ls("package:SeuratObject")) {
    seurat_obj <- JoinLayers(seurat_obj)
  }
}, error = function(e) { })

DefaultAssay(seurat_obj) <- "SCT"

# Key markers from original Sanghoon analysis
markers <- list(
  # Epithelial - Luminal/Cancer
  Epithelial_Luminal = c("Epcam", "Krt18", "Krt8"),
  # Epithelial - Myoepithelial/Basal
  Myoepithelial = c("Tp63", "Krt5", "Krt17", "Acta2", "Mylk", "Myh11"),
  # Endothelial
  Endothelial = c("Pecam1", "Cdh5", "Flt1", "Eng"),
  # Fibroblast
  Fibroblast = c("Col1a1", "Col1a2", "Pdgfrb", "Dcn"),
  # Myeloid/Macrophage
  Myeloid = c("Cd68", "Cd14", "Csf1r", "Adgre1", "Tyrobp"),
  # T/NK cells
  NKT = c("Cd3d", "Cd3e", "Cd3g", "Ncam1", "Klrd1"),
  # Pan-immune
  Immune = c("Ptprc")
)

# Check which markers exist
cat("Checking marker availability:\n")
all_markers <- unique(unlist(markers))
available <- all_markers[all_markers %in% rownames(seurat_obj)]
missing <- all_markers[!all_markers %in% rownames(seurat_obj)]
cat("  Available:", length(available), "/", length(all_markers), "\n")
if (length(missing) > 0) {
  cat("  Missing:", paste(missing, collapse = ", "), "\n")
}
cat("\n")

# Get expression data
expr_data <- GetAssayData(seurat_obj, layer = "data")

# Function to calculate mean expression by cell type
calc_mean_expr <- function(genes, cell_types) {
  genes <- genes[genes %in% rownames(expr_data)]
  if (length(genes) == 0) return(NULL)

  result <- sapply(unique(cell_types), function(ct) {
    cells <- names(cell_types)[cell_types == ct]
    mean(rowMeans(as.matrix(expr_data[genes, cells, drop = FALSE])))
  })
  return(result)
}

# Calculate expression by scType cell type
cell_types <- seurat_obj$sctype_cluster_celltype

cat("=== Mean Marker Expression by scType Cell Type ===\n\n")

expr_summary <- data.frame(CellType = unique(cell_types))

for (marker_group in names(markers)) {
  genes <- markers[[marker_group]]
  genes <- genes[genes %in% rownames(expr_data)]

  if (length(genes) > 0) {
    mean_expr <- sapply(unique(cell_types), function(ct) {
      cells <- names(cell_types)[cell_types == ct]
      mean(rowMeans(as.matrix(expr_data[genes, cells, drop = FALSE])))
    })
    expr_summary[[marker_group]] <- mean_expr[expr_summary$CellType]
  }
}

# Round for display
expr_summary[, -1] <- round(expr_summary[, -1], 3)

# Sort by a meaningful order
type_order <- c("CancerEpithelial", "Myoepithelial", "Endothelial", "Fibroblast",
                "Monocyte", "Macrophage", "M1_Macrophage", "M2_Macrophage",
                "DendriticCell", "NKcell", "CD4Tcell", "CD8Tcell", "Treg", "NaiveTcell")
expr_summary <- expr_summary %>%
  mutate(CellType = factor(CellType, levels = type_order)) %>%
  arrange(CellType)

print(expr_summary)

cat("\n=== Spot Check Analysis ===\n\n")

# 1. Endothelial validation
cat("1. ENDOTHELIAL VALIDATION (26.7% of cells)\n")
cat("   Expected: HIGH Pecam1, Cdh5, Flt1; LOW Epcam, Krt18\n")

endo_cells <- names(cell_types)[cell_types == "Endothelial"]
if (length(endo_cells) > 0) {
  endo_markers <- c("Pecam1", "Cdh5", "Flt1", "Epcam", "Krt18", "Cd68", "Ptprc")
  endo_markers <- endo_markers[endo_markers %in% rownames(expr_data)]

  endo_expr <- rowMeans(as.matrix(expr_data[endo_markers, endo_cells]))
  cat("   Mean expression in 'Endothelial' cells:\n")
  print(round(sort(endo_expr, decreasing = TRUE), 3))

  # Compare to CancerEpithelial
  epi_cells <- names(cell_types)[cell_types == "CancerEpithelial"]
  if (length(epi_cells) > 0) {
    epi_expr <- rowMeans(as.matrix(expr_data[endo_markers, epi_cells]))
    cat("\n   Comparison - CancerEpithelial cells:\n")
    print(round(sort(epi_expr, decreasing = TRUE), 3))
  }
}

cat("\n2. MYOEPITHELIAL VALIDATION (16.7% of cells)\n")
cat("   Expected: HIGH Tp63, Krt5, Acta2; MODERATE Krt18\n")

myo_cells <- names(cell_types)[cell_types == "Myoepithelial"]
if (length(myo_cells) > 0) {
  myo_markers <- c("Tp63", "Krt5", "Krt17", "Acta2", "Mylk", "Epcam", "Krt18")
  myo_markers <- myo_markers[myo_markers %in% rownames(expr_data)]

  myo_expr <- rowMeans(as.matrix(expr_data[myo_markers, myo_cells]))
  cat("   Mean expression in 'Myoepithelial' cells:\n")
  print(round(sort(myo_expr, decreasing = TRUE), 3))
}

cat("\n3. IMMUNE CELL VALIDATION (8.67% of cells)\n")
cat("   Expected: HIGH Ptprc (CD45), Cd68/Cd14 (myeloid), Cd3d (T cells)\n")

# Myeloid
mono_cells <- names(cell_types)[cell_types == "Monocyte"]
if (length(mono_cells) > 0) {
  immune_markers <- c("Ptprc", "Cd68", "Cd14", "Csf1r", "Cd3d", "Epcam", "Krt18")
  immune_markers <- immune_markers[immune_markers %in% rownames(expr_data)]

  mono_expr <- rowMeans(as.matrix(expr_data[immune_markers, mono_cells]))
  cat("   Monocyte cells:\n")
  print(round(sort(mono_expr, decreasing = TRUE), 3))
}

# T cells
cd4_cells <- names(cell_types)[cell_types == "CD4Tcell"]
if (length(cd4_cells) > 0) {
  t_markers <- c("Ptprc", "Cd3d", "Cd3e", "Stat4", "Cd68", "Epcam")
  t_markers <- t_markers[t_markers %in% rownames(expr_data)]

  cd4_expr <- rowMeans(as.matrix(expr_data[t_markers, cd4_cells]))
  cat("\n   CD4Tcell cells:\n")
  print(round(sort(cd4_expr, decreasing = TRUE), 3))
}

cat("\n4. CLUSTER-LEVEL SPOT CHECK\n")
cat("   Checking specific clusters with discrepancies...\n\n")

# Check clusters that scType called Endothelial but original called CancerEpithelial
discrepant_clusters <- c("1", "2", "8", "14", "15", "16")

for (cl in discrepant_clusters) {
  cl_cells <- colnames(seurat_obj)[seurat_obj$seurat_clusters == cl]
  if (length(cl_cells) == 0) next

  sctype_call <- unique(seurat_obj$sctype_cluster_celltype[seurat_obj$seurat_clusters == cl])

  key_markers <- c("Epcam", "Krt18", "Pecam1", "Cdh5", "Ptprc", "Cd68", "Col1a1")
  key_markers <- key_markers[key_markers %in% rownames(expr_data)]

  cl_expr <- rowMeans(as.matrix(expr_data[key_markers, cl_cells]))

  cat(sprintf("   Cluster %s (n=%d, scType=%s):\n", cl, length(cl_cells), sctype_call))
  print(round(sort(cl_expr, decreasing = TRUE), 3))
  cat("\n")
}

cat("=== Interpretation ===\n\n")
cat("If Endothelial cells show:\n")
cat("  - HIGH Pecam1/Cdh5/Flt1 and LOW Epcam/Krt18 → Correctly called\n")
cat("  - HIGH Epcam/Krt18 → Possibly misclassified epithelial\n\n")

cat("If Myoepithelial cells show:\n")
cat("  - HIGH Tp63/Krt5/Acta2 → Correctly called\n")
cat("  - Only HIGH Krt18, LOW Tp63 → Possibly misclassified luminal\n\n")

cat("=== Done ===\n")
