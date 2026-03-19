#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/03b_sctype_annotate.R
# Cell-by-cell annotation using scType
#
# This script provides more accurate cell type annotation by scoring
# individual cells rather than cluster averages, which prevents rare
# cell types (like immune cells) from being drowned out by dominant
# epithelial populations.
#
# Inputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_integrated.rds
#
# Outputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_annotated_sctype.rds
#   - analysis/02_rat_snrnaseq/outputs/sctype_scores.csv
#   - figures/by_analysis/rat_snrnaseq/sctype_annotation_umap.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(HGNChelper)  # For gene symbol checking
  library(openxlsx)    # For reading Excel marker database
})

# -----------------------------------------------------------------------------
# scType functions (from https://github.com/IanevskiAleksandr/sc-type)
# -----------------------------------------------------------------------------

# Gene sets preparation function
gene_sets_prepare <- function(path_to_db_file, cell_type) {
  # Read marker database
  if (grepl("^http", path_to_db_file)) {
    cell_markers <- openxlsx::read.xlsx(path_to_db_file)
  } else {
    cell_markers <- openxlsx::read.xlsx(path_to_db_file)
  }

  cell_markers <- cell_markers[cell_markers$tissueType == cell_type, ]
  cell_markers$geneSymbolmore1 <- gsub(" ", "", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 <- gsub(" ", "", cell_markers$geneSymbolmore2)

  # Prepare positive markers
  gs_positive <- lapply(1:nrow(cell_markers), function(j) {
    markers <- gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[j], ",")))
    markers <- markers[markers != "NA" & markers != ""]
    markers
  })
  names(gs_positive) <- cell_markers$cellName

  # Prepare negative markers
  gs_negative <- lapply(1:nrow(cell_markers), function(j) {
    markers <- gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[j], ",")))
    markers <- markers[markers != "NA" & markers != ""]
    markers
  })
  names(gs_negative) <- cell_markers$cellName

  return(list(gs_positive = gs_positive, gs_negative = gs_negative))
}

# scType scoring function
sctype_score <- function(scRNAseqData, scaled = TRUE, gs, gs2 = NULL,
                          gene_names_to_uppercase = FALSE, ...) {

  # Check input
  if (is.null(gs)) {
    stop("Please provide gene sets (gs)")
  }

  # Gene names to uppercase
  if (gene_names_to_uppercase) {
    rownames(scRNAseqData) <- toupper(rownames(scRNAseqData))
  }

  # Subselect genes only found in data
  names_gs_cp <- names(gs)
  names_gs_2_cp <- names(gs2)

  gs <- lapply(1:length(gs), function(d_) {
    GeneIn <- rownames(scRNAseqData)[rownames(scRNAseqData) %in% gs[[d_]]]
    GeneIn
  })
  names(gs) <- names_gs_cp

  gs2 <- lapply(1:length(gs2), function(d_) {
    GeneIn <- rownames(scRNAseqData)[rownames(scRNAseqData) %in% gs2[[d_]]]
    GeneIn
  })
  names(gs2) <- names_gs_2_cp

  # Compute marker specificity
  cell_markers_genes_score <- sort(unique(unlist(gs)))

  # Scale data if not already scaled
  if (!scaled) {
    Z <- t(scale(t(scRNAseqData)))
  } else {
    Z <- scRNAseqData
  }

  # Calculate scores for each cell type
  es <- sapply(1:length(gs), function(gss_) {
    gs_genes <- gs[[gss_]]
    gs2_genes <- gs2[[gss_]]

    # Calculate positive score
    if (length(gs_genes) > 0) {
      es_gene <- Z[gs_genes, , drop = FALSE]
      sum_pos <- colSums(es_gene) / sqrt(length(gs_genes))
    } else {
      sum_pos <- 0
    }

    # Calculate negative score
    if (!is.null(gs2_genes) && length(gs2_genes) > 0) {
      es_gene2 <- Z[gs2_genes, , drop = FALSE]
      sum_neg <- colSums(es_gene2) / sqrt(length(gs2_genes))
    } else {
      sum_neg <- 0
    }

    # Combined score
    sum_pos - sum_neg
  })

  if (is.vector(es)) {
    es <- matrix(es, nrow = 1)
  }

  colnames(es) <- names(gs)
  rownames(es) <- colnames(scRNAseqData)

  return(t(es))
}

# -----------------------------------------------------------------------------
# Custom Rat Mammary Marker Database
# -----------------------------------------------------------------------------
# Based on the original Sanghoon analysis markers and Li et al. 2020 Cell Reports

create_rat_mammary_markers <- function() {
  markers <- data.frame(
    tissueType = rep("Rat_Mammary", 14),
    cellName = c(
      "CancerEpithelial",
      "Myoepithelial",
      "Fibroblast",
      "Endothelial",
      "Monocyte",
      "Macrophage",
      "M1_Macrophage",
      "M2_Macrophage",
      "DendriticCell",
      "NKcell",
      "CD4Tcell",
      "CD8Tcell",
      "Treg",
      "NaiveTcell"
    ),
    # Positive markers (genes expressed in this cell type)
    geneSymbolmore1 = c(
      # CancerEpithelial - luminal markers
      "Epcam,Krt18,Krt8,Cd24,Prlr,Cited1,Pgr,Esr1",
      # Myoepithelial
      "Tp63,Krt17,Krt5,Acta2,Mylk,Myh11",
      # Fibroblast
      "Col1a1,Col1a2,Col3a1,Fn1,Pdgfrb,Dcn,C1r,Fap",
      # Endothelial
      "Pecam1,Cdh5,Flt1,Eng,Sox17,Sele,Ramp2",
      # Monocyte
      "Cd14,Csf1r,Adgre1,Ace,Fcgr3a,Cd68,Tyrobp",
      # Macrophage
      "Cd68,Apoe,Fabp5,Egr1,Cx3cr1,Slc2a1,Csf1r",
      # M1_Macrophage
      "Cd68,Cd86,Cd80,Nos2,Il1b,Tnf",
      # M2_Macrophage
      "Cd68,Msr1,Mrc1,Cd163,Clec10a,Arg1",
      # DendriticCell
      "Flt3,Gzmb,Itgax,Thbd,Fscn1,Cd83,Irf7,Lamp3,Clec9a",
      # NKcell
      "Ncam1,Klrd1,Il2rb,Klrb1,Klrc1,Klrk1,Ncr1,Areg,Xcl1,Gzmb",
      # CD4Tcell
      "Cd3d,Cd3e,Cd3g,Stat4,Il7r,Cd40lg,Btla",
      # CD8Tcell
      "Cd3d,Cd3e,Cd3g,Cd8a,Cd8b,Gzma,Gzmb",
      # Treg
      "Cd3d,Cd3e,Foxp3,Il2ra,Ctla4",
      # NaiveTcell
      "Cd3d,Cd3e,Ccr7,Sell,Il7r,Cd27"
    ),
    # Negative markers (genes NOT expressed in this cell type)
    geneSymbolmore2 = c(
      # CancerEpithelial - NOT immune markers
      "Ptprc,Cd3d,Cd68,Pecam1,Col1a1",
      # Myoepithelial - NOT luminal
      "Esr1,Pgr,Prlr,Cited1",
      # Fibroblast - NOT epithelial or immune
      "Epcam,Krt18,Ptprc,Cd68,Pecam1",
      # Endothelial - NOT epithelial or immune
      "Epcam,Krt18,Ptprc,Cd68,Col1a1",
      # Monocyte - NOT epithelial, NOT T cell
      "Epcam,Krt18,Cd3d,Cd3e,Ncam1",
      # Macrophage - NOT epithelial, NOT T cell
      "Epcam,Krt18,Cd3d,Cd3e,Ncam1,Pecam1",
      # M1_Macrophage - NOT M2 markers
      "Mrc1,Msr1,Cd163,Arg1",
      # M2_Macrophage - NOT M1 markers
      "Cd86,Nos2,Il1b",
      # DendriticCell - NOT T cell or macrophage
      "Cd3d,Cd3e,Mrc1,Msr1",
      # NKcell - NOT T cell markers
      "Cd3d,Cd3e,Cd4,Cd8a,Cd8b,Foxp3",
      # CD4Tcell - NOT CD8 or NK
      "Cd8a,Cd8b,Ncam1,Klrd1,Foxp3",
      # CD8Tcell - NOT CD4
      "Stat4,Foxp3,Ncam1",
      # Treg - NOT CD8 or NK
      "Cd8a,Cd8b,Ncam1,Klrd1",
      # NaiveTcell - NOT effector
      "Gzma,Gzmb,Foxp3"
    ),
    stringsAsFactors = FALSE
  )

  return(markers)
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

check_file_exists <- function(filepath, description = "file") {
  if (!file.exists(filepath)) {
    stop(sprintf("ERROR: %s not found: %s", description, filepath))
  }
  cat(sprintf("  Found: %s\n", basename(filepath)))
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}

# -----------------------------------------------------------------------------
# Main Script
# -----------------------------------------------------------------------------

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
fig_dir <- file.path(project_root, "figures/by_analysis/rat_snrnaseq")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Rat snRNA-seq scType Cell Annotation ===\n")
cat("Random seed: 12345\n")
cat("Project root:", project_root, "\n")
cat("Output directory:", output_dir, "\n\n")

# Input file - use integrated object (before previous annotation)
input_rds <- file.path(output_dir, "seurat_integrated.rds")

cat("Checking input files...\n")
check_file_exists(input_rds, "Integrated Seurat object")
cat("\n")

# -----------------------------------------------------------------------------
# Step 1: Load Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading integrated Seurat object...\n")
seurat_obj <- readRDS(input_rds)

if (!inherits(seurat_obj, "Seurat")) {
  stop("ERROR: Loaded object is not a valid Seurat object")
}

cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Features:", nrow(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Run UMAP and Clustering (match original parameters)
# -----------------------------------------------------------------------------
cat("\nStep 2: Running UMAP and clustering...\n")

# Match original analysis: resolution 1.5, 40 clusters
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30,
                       seed.use = 12345, verbose = FALSE)

seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30,
                             k.param = 15, verbose = FALSE)

# Use resolution 1.5 to match original analysis (40 clusters)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.5,
                            random.seed = 12345, verbose = FALSE)

n_clusters <- length(unique(Idents(seurat_obj)))
cat("  UMAP complete\n")
cat("  Found", n_clusters, "clusters (resolution 1.5)\n")

# -----------------------------------------------------------------------------
# Step 3: Prepare scType marker database
# -----------------------------------------------------------------------------
cat("\nStep 3: Preparing rat mammary marker database...\n")

# Create custom marker database
markers_db <- create_rat_mammary_markers()
cat("  Defined", nrow(markers_db), "cell types\n")

# Save marker database for reference
write.csv(markers_db, file.path(output_dir, "sctype_marker_database.csv"), row.names = FALSE)

# Prepare gene sets
gs_list <- list(gs_positive = list(), gs_negative = list())

for (i in 1:nrow(markers_db)) {
  cell_type <- markers_db$cellName[i]

  # Positive markers
  pos_markers <- unlist(strsplit(markers_db$geneSymbolmore1[i], ","))
  pos_markers <- trimws(pos_markers)
  pos_markers <- pos_markers[pos_markers != "" & pos_markers != "NA"]
  gs_list$gs_positive[[cell_type]] <- pos_markers

  # Negative markers
  neg_markers <- unlist(strsplit(markers_db$geneSymbolmore2[i], ","))
  neg_markers <- trimws(neg_markers)
  neg_markers <- neg_markers[neg_markers != "" & neg_markers != "NA"]
  gs_list$gs_negative[[cell_type]] <- neg_markers
}

# Check marker coverage
all_pos_markers <- unique(unlist(gs_list$gs_positive))
all_neg_markers <- unique(unlist(gs_list$gs_negative))
genes_in_data <- rownames(seurat_obj)

pos_found <- sum(all_pos_markers %in% genes_in_data)
neg_found <- sum(all_neg_markers %in% genes_in_data)

cat("  Positive markers: ", pos_found, "/", length(all_pos_markers), " found in data\n", sep = "")
cat("  Negative markers: ", neg_found, "/", length(all_neg_markers), " found in data\n", sep = "")

# Show missing markers
missing_pos <- all_pos_markers[!all_pos_markers %in% genes_in_data]
if (length(missing_pos) > 0) {
  cat("  Missing positive markers:", paste(head(missing_pos, 10), collapse = ", "),
      ifelse(length(missing_pos) > 10, "...", ""), "\n")
}

# -----------------------------------------------------------------------------
# Step 4: Run scType scoring (cell-by-cell)
# -----------------------------------------------------------------------------
cat("\nStep 4: Running scType cell-by-cell scoring...\n")
# Seurat v5 may require joining layers, but v4-created objects may not support it
tryCatch({
  if ("JoinLayers" %in% ls("package:SeuratObject")) {
    seurat_obj <- JoinLayers(seurat_obj)
    cat("  Joined Seurat v5 layers\n")
  }
}, error = function(e) {
  cat("  Note: JoinLayers not needed for this Seurat object (likely v4 format)\n")
})


# Get scaled data
# First ensure data is scaled
DefaultAssay(seurat_obj) <- "SCT"

# Try to get scale.data, if not available, scale it
tryCatch({
  scaled_data <- GetAssayData(seurat_obj, layer = "scale.data")
  if (ncol(scaled_data) == 0 || nrow(scaled_data) == 0) {
    stop("Empty scale data")
  }
  cat("  Using existing scaled data\n")
}, error = function(e) {
  cat("  Scaling data...\n")
  seurat_obj <<- ScaleData(seurat_obj, verbose = FALSE)
  scaled_data <<- GetAssayData(seurat_obj, layer = "scale.data")
})

scaled_data <- as.matrix(GetAssayData(seurat_obj, layer = "scale.data"))
cat("  Scaled data matrix:", nrow(scaled_data), "genes x", ncol(scaled_data), "cells\n")

# Run scType scoring
cat("  Computing scType scores for each cell...\n")

# Compute scores
es.max <- matrix(0, nrow = length(gs_list$gs_positive), ncol = ncol(scaled_data))
rownames(es.max) <- names(gs_list$gs_positive)
colnames(es.max) <- colnames(scaled_data)

for (cell_type in names(gs_list$gs_positive)) {
  # Get markers found in data
  pos_genes <- gs_list$gs_positive[[cell_type]]
  pos_genes <- pos_genes[pos_genes %in% rownames(scaled_data)]

  neg_genes <- gs_list$gs_negative[[cell_type]]
  neg_genes <- neg_genes[neg_genes %in% rownames(scaled_data)]

  if (length(pos_genes) > 0) {
    # Positive score: sum of expression / sqrt(n_markers)
    pos_score <- colSums(scaled_data[pos_genes, , drop = FALSE]) / sqrt(length(pos_genes))
  } else {
    pos_score <- 0
  }

  if (length(neg_genes) > 0) {
    # Negative score: penalize if negative markers are expressed
    neg_score <- colSums(scaled_data[neg_genes, , drop = FALSE]) / sqrt(length(neg_genes))
  } else {
    neg_score <- 0
  }

  # Final score = positive - negative
  es.max[cell_type, ] <- pos_score - neg_score
}

cat("  Scoring complete\n")

# -----------------------------------------------------------------------------
# Step 5: Assign cell types
# -----------------------------------------------------------------------------
cat("\nStep 5: Assigning cell types to cells...\n")

# For each cell, find the cell type with maximum score
cell_assignments <- apply(es.max, 2, function(x) {
  if (max(x) > 0) {
    names(x)[which.max(x)]
  } else {
    "Unknown"
  }
})

# Create a confidence score (difference between top 2 scores)
confidence_scores <- apply(es.max, 2, function(x) {
  sorted <- sort(x, decreasing = TRUE)
  if (length(sorted) >= 2) {
    sorted[1] - sorted[2]
  } else {
    sorted[1]
  }
})
# Name the assignments with cell barcodes
names(cell_assignments) <- Cells(seurat_obj)
names(confidence_scores) <- Cells(seurat_obj)


# Add to Seurat object
seurat_obj$sctype_celltype <- cell_assignments
seurat_obj$sctype_confidence <- confidence_scores

# Summary
cat("\nCell type distribution:\n")
cell_counts <- table(seurat_obj$sctype_celltype)
print(sort(cell_counts, decreasing = TRUE))

# Calculate percentages
cat("\nPercentages:\n")
cell_pct <- round(100 * cell_counts / sum(cell_counts), 2)
print(sort(cell_pct, decreasing = TRUE))

# -----------------------------------------------------------------------------
# Step 6: Create simplified main cell type categories
# -----------------------------------------------------------------------------
cat("\nStep 6: Creating main cell type categories...\n")

# Map detailed types to main categories (matching original analysis naming)
main_type_mapping <- c(
  "CancerEpithelial" = "CancerEpithelial",
  "Myoepithelial" = "Myoepithelial",
  "Fibroblast" = "Fibroblast",
  "Endothelial" = "Endothelial",
  "Monocyte" = "Myeloid",
  "Macrophage" = "Myeloid",
  "M1_Macrophage" = "Myeloid",
  "M2_Macrophage" = "Myeloid",
  "DendriticCell" = "DendriticCell",
  "NKcell" = "NKTcell",
  "CD4Tcell" = "NKTcell",
  "CD8Tcell" = "NKTcell",
  "Treg" = "NKTcell",
  "NaiveTcell" = "NKTcell",
  "Unknown" = "Unknown"
)

# Create mapped cell types with proper naming
mapped_types <- main_type_mapping[seurat_obj$sctype_celltype]
names(mapped_types) <- Cells(seurat_obj)
seurat_obj$CellTypeByMarker_RatsnRNAseq <- mapped_types
seurat_obj$CellTypeMacroTcell_RatsnRNAseq <- seurat_obj$sctype_celltype

cat("\nMain cell type distribution:\n")
main_counts <- table(seurat_obj$CellTypeByMarker_RatsnRNAseq)
print(sort(main_counts, decreasing = TRUE))

# -----------------------------------------------------------------------------
# Step 7: Generate Visualizations
# -----------------------------------------------------------------------------
cat("\nStep 7: Generating visualizations...\n")

# UMAP by scType cell type
pdf(file.path(fig_dir, "sctype_annotation_umap.pdf"), width = 14, height = 10)

# Detailed cell types
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "sctype_celltype",
              label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("scType Cell-by-Cell Annotation (Detailed)") +
  theme(legend.position = "right")
print(p1)

# Main cell types
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "CellTypeByMarker_RatsnRNAseq",
              label = TRUE, label.size = 5, repel = TRUE) +
  ggtitle("scType Cell-by-Cell Annotation (Main Types)") +
  theme(legend.position = "right")
print(p2)

# By cluster
p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, label.size = 4) +
  ggtitle(paste0("Seurat Clusters (n=", n_clusters, ", resolution 1.5)"))
print(p3)

# By age group
p4 <- DimPlot(seurat_obj, reduction = "umap", group.by = "AgeGroup",
              label = FALSE) +
  ggtitle("Age Group Distribution")
print(p4)

# Confidence score
p5 <- FeaturePlot(seurat_obj, features = "sctype_confidence",
                   cols = c("lightgrey", "blue")) +
  ggtitle("scType Confidence Score")
print(p5)

dev.off()

# Save individual plots as PNG and SVG
ggsave(file.path(fig_dir, "sctype_detailed_umap.png"), p1, width = 14, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "sctype_detailed_umap.svg"), p1, width = 14, height = 10)
ggsave(file.path(fig_dir, "sctype_main_umap.png"), p2, width = 12, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "sctype_main_umap.svg"), p2, width = 12, height = 10)

cat("  Plots saved to:", fig_dir, "\n")

# -----------------------------------------------------------------------------
# Step 8: Compare with cluster-level assignments
# -----------------------------------------------------------------------------
cat("\nStep 8: Cluster-level cell type composition...\n")

# For each cluster, what's the dominant cell type?
cluster_composition <- seurat_obj@meta.data %>%
  group_by(seurat_clusters, sctype_celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(
    total = sum(n),
    pct = round(100 * n / total, 1)
  ) %>%
  arrange(seurat_clusters, desc(pct))

# Dominant type per cluster
dominant_per_cluster <- cluster_composition %>%
  group_by(seurat_clusters) %>%
  slice_max(pct, n = 1) %>%
  select(seurat_clusters, dominant_type = sctype_celltype, pct)

cat("\nDominant cell type per cluster:\n")
print(as.data.frame(dominant_per_cluster), row.names = FALSE)

# -----------------------------------------------------------------------------
# Step 9: Save Results
# -----------------------------------------------------------------------------
cat("\nStep 9: Saving results...\n")

# Save annotated Seurat object
saveRDS(seurat_obj, file.path(output_dir, "seurat_annotated_sctype.rds"))

# Save scType scores matrix
write.csv(t(es.max), file.path(output_dir, "sctype_scores.csv"))

# Save cell type assignments
assignments_df <- data.frame(
  cell_id = colnames(seurat_obj),
  cluster = seurat_obj$seurat_clusters,
  sctype_celltype = seurat_obj$sctype_celltype,
  main_celltype = seurat_obj$CellTypeByMarker_RatsnRNAseq,
  confidence = seurat_obj$sctype_confidence,
  AgeGroup = seurat_obj$AgeGroup,
  orig.ident = seurat_obj$orig.ident
)
write.csv(assignments_df, file.path(output_dir, "sctype_cell_assignments.csv"), row.names = FALSE)

# Save cluster composition
write.csv(cluster_composition, file.path(output_dir, "sctype_cluster_composition.csv"), row.names = FALSE)

cat("  Saved:\n")
cat("    - seurat_annotated_sctype.rds\n")
cat("    - sctype_scores.csv\n")
cat("    - sctype_cell_assignments.csv\n")
cat("    - sctype_cluster_composition.csv\n")
cat("    - sctype_marker_database.csv\n")

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Total cells:", ncol(seurat_obj), "\n")
cat("Clusters:", n_clusters, "\n")
cat("\nCell type distribution:\n")

summary_df <- seurat_obj@meta.data %>%
  group_by(CellTypeByMarker_RatsnRNAseq) %>%
  summarise(
    n = n(),
    pct = round(100 * n() / ncol(seurat_obj), 2)
  ) %>%
  arrange(desc(n))

print(as.data.frame(summary_df), row.names = FALSE)

# Immune cell summary
immune_types <- c("Myeloid", "DendriticCell", "NKTcell")
immune_n <- sum(seurat_obj$CellTypeByMarker_RatsnRNAseq %in% immune_types)
immune_pct <- round(100 * immune_n / ncol(seurat_obj), 2)

cat("\nImmune cells total:", immune_n, "(", immune_pct, "%)\n")

cat("\n=== scType annotation complete ===\n")
