#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/02_run_gsva.R
# Run Gene Set Variation Analysis on estrogen-related pathways
#
# Inputs:
#   - analysis/01_human_bulk_rnaseq/outputs/vst_normalized_matrix.rds
#   - data/human_bulk_rnaseq/external/SuppleTable1_UpregulatedByE1_NotE2_407g.txt
#
# Outputs:
#   - analysis/01_human_bulk_rnaseq/outputs/gsva_estrogen_pathways.rds
#   - analysis/01_human_bulk_rnaseq/outputs/gsva_heatmap.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(GSVA)
  library(msigdbr)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# Helper function to check file existence
check_file_exists <- function(filepath, description = "file") {
  if (!file.exists(filepath)) {
    stop(sprintf("ERROR: %s not found: %s", description, filepath))
  }
  cat(sprintf("  Found: %s\n", basename(filepath)))
}

# Define paths
script_dir <- dirname(sys.frame(1)$ofile)
project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
data_dir <- file.path(project_root, "data/human_bulk_rnaseq")

cat("=== GSVA Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load Normalized Data
# -----------------------------------------------------------------------------
cat("Step 1: Loading normalized data...\n")
vst_file <- file.path(output_dir, "vst_normalized_matrix.rds")
check_file_exists(vst_file, "VST normalized matrix")
vst_matrix <- readRDS(vst_file)
vst_mat <- vst_matrix %>%
  tibble::column_to_rownames("GeneSymb") %>%
  as.matrix()
cat("  Matrix:", nrow(vst_mat), "genes x", ncol(vst_mat), "samples\n")

# -----------------------------------------------------------------------------
# Step 2: Load Gene Sets
# -----------------------------------------------------------------------------
cat("Step 2: Loading gene sets...\n")

# E1-upregulated genes
e1_file <- file.path(data_dir, "external/SuppleTable1_UpregulatedByE1_NotE2_407g.txt")
if (file.exists(e1_file)) {
  e1_genes <- fread(e1_file, header = TRUE, stringsAsFactors = FALSE)
  colnames(e1_genes) <- gsub(" ", "", colnames(e1_genes))
  e1_genes$Gene <- gsub("-.*", "", e1_genes$Gene)
  e1_genes <- e1_genes %>% filter(!duplicated(Gene))
  e1_gene_list <- list(E1UpRegGene = e1_genes$Gene)
  cat("  E1-upregulated genes:", length(e1_gene_list$E1UpRegGene), "\n")
} else {
  e1_gene_list <- list()
  cat("  Warning: E1 gene file not found\n")
}

# Hallmark estrogen pathways
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
hallmark_estrogen <- hallmark_list[c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"
)]

# Reactome estrogen pathway
reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
reactome_list <- split(reactome_sets$gene_symbol, reactome_sets$gs_name)
reactome_estrogen <- reactome_list["REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION"]

# WikiPathways
wiki_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "WIKIPATHWAYS")
wiki_list <- split(wiki_sets$gene_symbol, wiki_sets$gs_name)
wiki_estrogen <- wiki_list["WP_ESTROGEN_SIGNALING_PATHWAY"]

# GO BP estrogen pathways
gobp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gobp_list <- split(gobp_sets$gene_symbol, gobp_sets$gs_name)
gobp_estrogen <- gobp_list[c(
  "GOBP_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS"
)]

# Combine all estrogen pathways
estrogen_pathways <- c(hallmark_estrogen, e1_gene_list, reactome_estrogen, wiki_estrogen, gobp_estrogen)
estrogen_pathways <- estrogen_pathways[!sapply(estrogen_pathways, is.null)]
cat("  Total pathways:", length(estrogen_pathways), "\n")

# Log which gene sets were found
cat("  Gene sets included in GSVA:\n")
for (gs_name in names(estrogen_pathways)) {
  n_genes <- length(estrogen_pathways[[gs_name]])
  n_in_data <- sum(estrogen_pathways[[gs_name]] %in% rownames(vst_mat))
  cat(sprintf("    - %s: %d genes (%d found in data)\n", gs_name, n_genes, n_in_data))
}

# -----------------------------------------------------------------------------
# Step 3: Run GSVA
# -----------------------------------------------------------------------------
cat("Step 3: Running GSVA...\n")

gsva_result <- gsva(
  gsvaParam(
    vst_mat,
    estrogen_pathways,
    kcdf = "Gaussian",  # Appropriate for continuous (VST) data
    maxDiff = TRUE
  )
)

cat("  GSVA result:", nrow(gsva_result), "pathways x", ncol(gsva_result), "samples\n")

# -----------------------------------------------------------------------------
# Step 4: Generate Heatmap
# -----------------------------------------------------------------------------
cat("Step 4: Generating heatmap...\n")

# Extract sample metadata from column names
sample_info <- data.frame(
  Sample = colnames(gsva_result),
  AgeRange = gsub(".*_(.*)Tumor.*", "\\1", colnames(gsva_result)),
  Group = ifelse(grepl("TumorAdj", colnames(gsva_result)), "TumorAdj", "Tumor")
)
rownames(sample_info) <- sample_info$Sample

# Annotation colors
ann_colors <- list(
  AgeRange = c(Young = "#4DAF4A", Middle = "#377EB8", Elderly = "#E41A1C"),
  Group = c(Tumor = "#984EA3", TumorAdj = "#FF7F00")
)

pdf(file.path(output_dir, "gsva_heatmap.pdf"), width = 14, height = 8)
pheatmap(
  gsva_result,
  annotation_col = sample_info[, c("AgeRange", "Group")],
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  main = "GSVA Estrogen Pathway Scores",
  scale = "row"
)
dev.off()

# -----------------------------------------------------------------------------
# Step 5: Save Outputs
# -----------------------------------------------------------------------------
cat("Step 5: Saving outputs...\n")
saveRDS(gsva_result, file.path(output_dir, "gsva_estrogen_pathways.rds"))
saveRDS(estrogen_pathways, file.path(output_dir, "estrogen_pathway_genesets.rds"))

cat("\n=== GSVA complete ===\n")
