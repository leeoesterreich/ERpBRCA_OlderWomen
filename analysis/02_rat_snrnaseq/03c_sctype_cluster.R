#!/usr/bin/env Rscript
# analysis/02_rat_snrnaseq/03c_sctype_cluster.R
# scType annotation at CLUSTER level (not cell-by-cell)
#
# This aggregates cell-level scType scores by cluster, then assigns
# one cell type per cluster - similar to the original Sanghoon analysis
# which manually assigned cell types to each of 40 clusters.
#
# Inputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_integrated.rds
#
# Outputs:
#   - analysis/02_rat_snrnaseq/outputs/seurat_annotated_sctype_cluster.rds
#   - figures/by_analysis/rat_snrnaseq/sctype_cluster_annotation.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(HGNChelper)
  library(openxlsx)
})

# Set Arial as default font for all plots
library(showtext)
font_add("Arial", "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Arial.ttf")
showtext_auto()
theme_set(theme_bw(base_size = 14, base_family = "Arial"))

# -----------------------------------------------------------------------------
# Custom Rat Mammary Marker Database (same as 03b)
# -----------------------------------------------------------------------------

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
    geneSymbolmore1 = c(
      "Epcam,Krt18,Krt8,Cd24,Prlr,Cited1,Pgr,Esr1",
      "Tp63,Krt17,Krt5,Acta2,Mylk,Myh11",
      "Col1a1,Col1a2,Col3a1,Fn1,Pdgfrb,Dcn,C1r,Fap",
      # Endothelial: removed Flt1 (also expressed in epithelial), emphasized Pecam1/Cdh5
      "Pecam1,Cdh5,Eng,Sox17,Sele,Ramp2,Vwf,Tie1",
      "Cd14,Csf1r,Adgre1,Ace,Fcgr3a,Cd68,Tyrobp",
      "Cd68,Apoe,Fabp5,Egr1,Cx3cr1,Slc2a1,Csf1r",
      "Cd68,Cd86,Cd80,Nos2,Il1b,Tnf",
      "Cd68,Msr1,Mrc1,Cd163,Clec10a,Arg1",
      "Flt3,Gzmb,Itgax,Thbd,Fscn1,Cd83,Irf7,Lamp3,Clec9a",
      "Ncam1,Klrd1,Il2rb,Klrb1,Klrc1,Klrk1,Ncr1,Areg,Xcl1,Gzmb",
      "Cd3d,Cd3e,Cd3g,Stat4,Il7r,Cd40lg,Btla",
      "Cd3d,Cd3e,Cd3g,Cd8a,Cd8b,Gzma,Gzmb",
      "Cd3d,Cd3e,Foxp3,Il2ra,Ctla4",
      "Cd3d,Cd3e,Ccr7,Sell,Il7r,Cd27"
    ),
    geneSymbolmore2 = c(
      "Ptprc,Cd3d,Cd68,Pecam1,Col1a1",
      "Esr1,Pgr,Prlr,Cited1",
      "Epcam,Krt18,Ptprc,Cd68,Pecam1",
      # Endothelial: strongly penalize epithelial markers
      "Epcam,Krt18,Krt8,Krt5,Ptprc,Cd68,Col1a1,Cd3d",
      "Epcam,Krt18,Cd3d,Cd3e,Ncam1",
      "Epcam,Krt18,Cd3d,Cd3e,Ncam1,Pecam1",
      "Mrc1,Msr1,Cd163,Arg1",
      "Cd86,Nos2,Il1b",
      "Cd3d,Cd3e,Mrc1,Msr1",
      "Cd3d,Cd3e,Cd4,Cd8a,Cd8b,Foxp3",
      "Cd8a,Cd8b,Ncam1,Klrd1,Foxp3",
      "Stat4,Foxp3,Ncam1",
      "Cd8a,Cd8b,Ncam1,Klrd1",
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

cat("=== scType Cluster-Level Annotation ===\n")
cat("Random seed: 12345\n")
cat("Project root:", project_root, "\n\n")

# Input file
input_rds <- file.path(output_dir, "seurat_integrated.rds")

cat("Checking input files...\n")
check_file_exists(input_rds, "Integrated Seurat object")
cat("\n")

# -----------------------------------------------------------------------------
# Step 1: Load Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading integrated Seurat object...\n")
seurat_obj <- readRDS(input_rds)
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Features:", nrow(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Run UMAP and Clustering (resolution 1.5 = 40 clusters like original)
# -----------------------------------------------------------------------------
cat("\nStep 2: Running UMAP and clustering (resolution 1.5)...\n")

seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30,
                       seed.use = 12345, verbose = FALSE)

seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30,
                             k.param = 15, verbose = FALSE)

seurat_obj <- FindClusters(seurat_obj, resolution = 1.5,
                            random.seed = 12345, verbose = FALSE)

n_clusters <- length(unique(Idents(seurat_obj)))
cat("  Found", n_clusters, "clusters\n")

# -----------------------------------------------------------------------------
# Step 3: Prepare marker database
# -----------------------------------------------------------------------------
cat("\nStep 3: Preparing marker database...\n")

markers_db <- create_rat_mammary_markers()

gs_list <- list(gs_positive = list(), gs_negative = list())
for (i in 1:nrow(markers_db)) {
  cell_type <- markers_db$cellName[i]

  pos_markers <- unlist(strsplit(markers_db$geneSymbolmore1[i], ","))
  pos_markers <- trimws(pos_markers)
  pos_markers <- pos_markers[pos_markers != "" & pos_markers != "NA"]
  gs_list$gs_positive[[cell_type]] <- pos_markers

  neg_markers <- unlist(strsplit(markers_db$geneSymbolmore2[i], ","))
  neg_markers <- trimws(neg_markers)
  neg_markers <- neg_markers[neg_markers != "" & neg_markers != "NA"]
  gs_list$gs_negative[[cell_type]] <- neg_markers
}

cat("  Defined", length(gs_list$gs_positive), "cell types\n")

# -----------------------------------------------------------------------------
# Step 4: Get scaled data and compute cell-level scores
# -----------------------------------------------------------------------------
cat("\nStep 4: Computing scType scores per cell...\n")

# Handle Seurat v5 layers
tryCatch({
  if ("JoinLayers" %in% ls("package:SeuratObject")) {
    seurat_obj <- JoinLayers(seurat_obj)
    cat("  Joined Seurat v5 layers\n")
  }
}, error = function(e) {
  cat("  Note: JoinLayers not needed\n")
})

DefaultAssay(seurat_obj) <- "SCT"

# Get or create scaled data
tryCatch({
  scaled_data <- GetAssayData(seurat_obj, layer = "scale.data")
  if (ncol(scaled_data) == 0 || nrow(scaled_data) == 0) {
    stop("Empty")
  }
  cat("  Using existing scaled data\n")
}, error = function(e) {
  cat("  Scaling data...\n")
  seurat_obj <<- ScaleData(seurat_obj, verbose = FALSE)
})

scaled_data <- as.matrix(GetAssayData(seurat_obj, layer = "scale.data"))
cat("  Matrix:", nrow(scaled_data), "genes x", ncol(scaled_data), "cells\n")

# Compute cell-level scores
es.max <- matrix(0, nrow = length(gs_list$gs_positive), ncol = ncol(scaled_data))
rownames(es.max) <- names(gs_list$gs_positive)
colnames(es.max) <- colnames(scaled_data)

for (cell_type in names(gs_list$gs_positive)) {
  pos_genes <- gs_list$gs_positive[[cell_type]]
  pos_genes <- pos_genes[pos_genes %in% rownames(scaled_data)]

  neg_genes <- gs_list$gs_negative[[cell_type]]
  neg_genes <- neg_genes[neg_genes %in% rownames(scaled_data)]

  if (length(pos_genes) > 0) {
    pos_score <- colSums(scaled_data[pos_genes, , drop = FALSE]) / sqrt(length(pos_genes))
  } else {
    pos_score <- 0
  }

  if (length(neg_genes) > 0) {
    neg_score <- colSums(scaled_data[neg_genes, , drop = FALSE]) / sqrt(length(neg_genes))
  } else {
    neg_score <- 0
  }

  es.max[cell_type, ] <- pos_score - neg_score
}

cat("  Cell-level scoring complete\n")

# -----------------------------------------------------------------------------
# Step 5: Aggregate scores by cluster and assign cluster-level cell types
# -----------------------------------------------------------------------------
cat("\nStep 5: Aggregating scores by cluster...\n")

# Get cluster assignments
clusters <- as.character(Idents(seurat_obj))
names(clusters) <- colnames(seurat_obj)

# For each cluster, sum the scores across cells and pick the max
cluster_scores <- data.frame(
  cluster = character(),
  cell_type = character(),
  score = numeric(),
  n_cells = integer(),
  stringsAsFactors = FALSE
)

cluster_assignments <- list()

for (cl in sort(unique(clusters))) {
  cells_in_cluster <- names(clusters)[clusters == cl]
  n_cells <- length(cells_in_cluster)

  # Sum scores across cells in this cluster
  cluster_score_sums <- rowSums(es.max[, cells_in_cluster, drop = FALSE])

  # Normalize by number of cells
  cluster_score_means <- cluster_score_sums / n_cells

  # Find the cell type with highest score
  best_type <- names(cluster_score_means)[which.max(cluster_score_means)]
  best_score <- max(cluster_score_means)

  cluster_assignments[[cl]] <- best_type

  # Store all scores for this cluster
  for (ct in names(cluster_score_means)) {
    cluster_scores <- rbind(cluster_scores, data.frame(
      cluster = cl,
      cell_type = ct,
      score = cluster_score_means[ct],
      n_cells = n_cells,
      stringsAsFactors = FALSE
    ))
  }
}

# Show cluster assignments
cluster_assignment_df <- data.frame(
  cluster = names(cluster_assignments),
  assigned_type = unlist(cluster_assignments),
  stringsAsFactors = FALSE
)

# Get cell counts per cluster
cluster_sizes <- table(clusters)
cluster_assignment_df$n_cells <- as.integer(cluster_sizes[cluster_assignment_df$cluster])

cat("\nCluster assignments:\n")
print(cluster_assignment_df %>% arrange(as.numeric(cluster)))

# -----------------------------------------------------------------------------
# Step 6: Assign cell types based on cluster membership
# -----------------------------------------------------------------------------
cat("\nStep 6: Assigning cell types to cells based on cluster...\n")

# Map cells to their cluster's assigned type
cell_types_by_cluster <- cluster_assignments[clusters]
names(cell_types_by_cluster) <- names(clusters)

seurat_obj$sctype_cluster_celltype <- unlist(cell_types_by_cluster)

# Create main type mapping
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
  "NaiveTcell" = "NKTcell"
)

mapped_main <- main_type_mapping[seurat_obj$sctype_cluster_celltype]
names(mapped_main) <- Cells(seurat_obj)
seurat_obj$CellTypeByMarker_RatsnRNAseq <- mapped_main
seurat_obj$CellTypeMacroTcell_RatsnRNAseq <- seurat_obj$sctype_cluster_celltype

# Summary
cat("\nDetailed cell type distribution:\n")
detail_counts <- table(seurat_obj$sctype_cluster_celltype)
print(sort(detail_counts, decreasing = TRUE))

cat("\nMain cell type distribution:\n")
main_counts <- table(seurat_obj$CellTypeByMarker_RatsnRNAseq)
print(sort(main_counts, decreasing = TRUE))

# Calculate percentages
cat("\nPercentages:\n")
main_pct <- round(100 * main_counts / sum(main_counts), 2)
print(sort(main_pct, decreasing = TRUE))

# Immune summary
immune_types <- c("Myeloid", "DendriticCell", "NKTcell")
immune_n <- sum(seurat_obj$CellTypeByMarker_RatsnRNAseq %in% immune_types)
immune_pct <- round(100 * immune_n / ncol(seurat_obj), 2)
cat("\nTotal immune cells:", immune_n, "(", immune_pct, "%)\n")

# -----------------------------------------------------------------------------
# Step 7: Compare with original Sanghoon cluster assignments
# -----------------------------------------------------------------------------
cat("\nStep 7: Comparison with original analysis...\n")

# Original Sanghoon assignments for 40 clusters (resolution 1.5)
original_assignments <- c(
  "0" = "CancerEpithelial", "1" = "CancerEpithelial", "2" = "CancerEpithelial",
  "3" = "CancerEpithelial", "4" = "CancerEpithelial", "5" = "CancerEpithelial",
  "6" = "CancerEpithelial", "7" = "CancerEpithelial", "8" = "CancerEpithelial",
  "9" = "Myoepithelial", "10" = "CancerEpithelial", "11" = "CancerEpithelial",
  "12" = "CancerEpithelial", "13" = "CancerEpithelial", "14" = "CancerEpithelial",
  "15" = "Myeloid", "16" = "NKTcell", "17" = "CancerEpithelial",
  "18" = "CancerEpithelial", "19" = "CancerEpithelial", "20" = "Myoepithelial",
  "21" = "Myeloid", "22" = "NKTcell", "23" = "CancerEpithelial",
  "24" = "CancerEpithelial", "25" = "Fibroblast", "26" = "CancerEpithelial",
  "27" = "DendriticCell", "28" = "Myeloid", "29" = "CancerEpithelial",
  "30" = "Endothelial", "31" = "CancerEpithelial", "32" = "Fibroblast",
  "33" = "Endothelial", "34" = "CancerEpithelial", "35" = "CancerEpithelial",
  "36" = "CancerEpithelial", "37" = "DendriticCell", "38" = "CancerEpithelial",
  "39" = "Myoepithelial"
)

# Compare our assignments to original
comparison_df <- cluster_assignment_df %>%
  mutate(
    original_type = ifelse(cluster %in% names(original_assignments),
                           original_assignments[cluster], NA),
    # Map our detailed types to main categories for comparison
    our_main_type = main_type_mapping[assigned_type],
    match = (our_main_type == original_type)
  )

cat("\nCluster-by-cluster comparison (scType vs Original):\n")
print(comparison_df %>%
        select(cluster, n_cells, assigned_type, our_main_type, original_type, match) %>%
        arrange(as.numeric(cluster)))

n_match <- sum(comparison_df$match, na.rm = TRUE)
n_total <- sum(!is.na(comparison_df$match))
cat("\nAgreement:", n_match, "/", n_total, "clusters match (",
    round(100 * n_match / n_total, 1), "%)\n")

# -----------------------------------------------------------------------------
# Step 8: Generate Visualizations
# -----------------------------------------------------------------------------
cat("\nStep 8: Generating visualizations...\n")

pdf(file.path(fig_dir, "sctype_cluster_annotation.pdf"), width = 14, height = 10)

# Detailed cell types
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "sctype_cluster_celltype",
              label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle("scType Cluster-Level Annotation (Detailed)") +
  theme(legend.position = "right")
print(p1)

# Main cell types
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "CellTypeByMarker_RatsnRNAseq",
              label = TRUE, label.size = 5, repel = TRUE) +
  ggtitle("scType Cluster-Level Annotation (Main Types)") +
  theme(legend.position = "right")
print(p2)

# By cluster number
p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, label.size = 3) +
  ggtitle(paste0("Seurat Clusters (n=", n_clusters, ")"))
print(p3)

# By age group
p4 <- DimPlot(seurat_obj, reduction = "umap", group.by = "AgeGroup") +
  ggtitle("Age Group Distribution")
print(p4)

dev.off()

# Save PNGs and SVGs
ggsave(file.path(fig_dir, "sctype_cluster_detailed.png"), p1, width = 14, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "sctype_cluster_detailed.svg"), p1, width = 14, height = 10)
ggsave(file.path(fig_dir, "sctype_cluster_main.png"), p2, width = 12, height = 10, dpi = 300)
ggsave(file.path(fig_dir, "sctype_cluster_main.svg"), p2, width = 12, height = 10)

cat("  Plots saved\n")

# -----------------------------------------------------------------------------
# Step 9: Save Results
# -----------------------------------------------------------------------------
cat("\nStep 9: Saving results...\n")

saveRDS(seurat_obj, file.path(output_dir, "seurat_annotated_sctype_cluster.rds"))
write.csv(cluster_assignment_df, file.path(output_dir, "sctype_cluster_assignments.csv"), row.names = FALSE)
write.csv(cluster_scores, file.path(output_dir, "sctype_cluster_scores.csv"), row.names = FALSE)
write.csv(comparison_df, file.path(output_dir, "sctype_vs_original_comparison.csv"), row.names = FALSE)

cat("  Saved:\n")
cat("    - seurat_annotated_sctype_cluster.rds\n")
cat("    - sctype_cluster_assignments.csv\n")
cat("    - sctype_cluster_scores.csv\n")
cat("    - sctype_vs_original_comparison.csv\n")

# -----------------------------------------------------------------------------
# Final Summary
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat("Total cells:", ncol(seurat_obj), "\n")
cat("Clusters:", n_clusters, "\n")
cat("Immune cells:", immune_n, "(", immune_pct, "%)\n")

cat("\nComparison with original Sanghoon analysis:\n")
cat("  - Original found ~9.8% immune cells\n")
cat("  - scType cluster-level found", immune_pct, "% immune cells\n")
cat("  - Cluster agreement:", round(100 * n_match / n_total, 1), "%\n")

cat("\n=== scType cluster-level annotation complete ===\n")
