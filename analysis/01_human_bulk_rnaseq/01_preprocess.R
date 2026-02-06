#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/01_preprocess.R
# Preprocess human bulk RNA-seq data: read counts, filter genes, normalize
#
# Inputs:
#   - data/human_bulk_rnaseq/raw/HumanERpAge_39404g168s_FeatureCount.txt
#   - data/human_bulk_rnaseq/raw/HumanERpAge_BulkRNAseq_SampleInformation.txt
#   - data/human_bulk_rnaseq/external/Homo_sapiens.gene_info.txt
#
# Outputs:
#   - analysis/01_human_bulk_rnaseq/outputs/dds_norm_AllAgeGroup.rds
#   - analysis/01_human_bulk_rnaseq/outputs/vst_normalized_matrix.rds
#   - analysis/01_human_bulk_rnaseq/outputs/qc_plots.pdf

# Set random seed for reproducibility
set.seed(12345)

# Load libraries
suppressPackageStartupMessages({
  library(tidyr)
  library(DESeq2)
  library(dplyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
})

# Helper function to check file existence
check_file_exists <- function(filepath, description = "file") {
  if (!file.exists(filepath)) {
    stop(sprintf("ERROR: %s not found: %s", description, filepath))
  }
  cat(sprintf("  Found: %s\n", basename(filepath)))
}

# Define paths relative to project root
project_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))
data_dir <- file.path(project_root, "data/human_bulk_rnaseq")
output_dir <- file.path(dirname(sys.frame(1)$ofile), "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input files
count_file <- file.path(data_dir, "raw/HumanERpAge_39404g168s_FeatureCount.txt")
sample_file <- file.path(data_dir, "raw/HumanERpAge_BulkRNAseq_SampleInformation.txt")
gene_annot_file <- file.path(data_dir, "external/Homo_sapiens.gene_info.txt")

cat("=== Human Bulk RNA-seq Preprocessing ===\n")
cat("Project root:", project_root, "\n")
cat("Output directory:", output_dir, "\n\n")

# Verify input files exist
cat("Checking input files...\n")
check_file_exists(count_file, "Count data file")
check_file_exists(sample_file, "Sample annotation file")
check_file_exists(gene_annot_file, "Gene annotation file")
cat("\n")

# -----------------------------------------------------------------------------
# Step 1: Read Gene Annotation
# -----------------------------------------------------------------------------
cat("Step 1: Reading gene annotation...\n")
gene_annot <- fread(gene_annot_file, header = TRUE, stringsAsFactors = FALSE)
n_genes_raw <- nrow(gene_annot)
gene_annot_filtered <- gene_annot %>%
  dplyr::select(Symbol, type_of_gene) %>%
  dplyr::filter(type_of_gene == "protein-coding", !grepl("^LOC\\d+", Symbol))
cat("  Total genes in annotation:", n_genes_raw, "\n")
cat("  Protein-coding genes (after filter):", nrow(gene_annot_filtered), "\n")
cat("  Genes removed:", n_genes_raw - nrow(gene_annot_filtered), "\n")

# -----------------------------------------------------------------------------
# Step 2: Read Count Data
# -----------------------------------------------------------------------------
cat("Step 2: Reading count data...\n")
count_data <- fread(count_file, header = TRUE, stringsAsFactors = FALSE)
colnames(count_data)[1] <- "GeneSymb"
colnames(count_data) <- gsub("_LEE(.*)", "", colnames(count_data))
cat("  Raw genes:", nrow(count_data), "\n")
cat("  Samples:", ncol(count_data) - 1, "\n")

# Filter to protein-coding genes
n_genes_before_filter <- nrow(count_data)
count_data_filtered <- count_data %>%
  dplyr::filter(GeneSymb %in% gene_annot_filtered$Symbol)
cat("  Protein-coding genes in data:", nrow(count_data_filtered), "\n")
cat("  Genes removed (non-protein-coding):", n_genes_before_filter - nrow(count_data_filtered), "\n")

# -----------------------------------------------------------------------------
# Step 3: Read Sample Annotation
# -----------------------------------------------------------------------------
cat("Step 3: Reading sample annotation...\n")
sample_annot <- fread(sample_file, stringsAsFactors = FALSE, header = TRUE)
colnames(sample_annot)[2] <- "SampleName"
sample_annot <- sample_annot %>%
  dplyr::mutate(
    AgeRange_Group = paste0(AgeRange, "_", Group),
    SampleNameGroup = paste0(SampleName, "_", AgeRange, Group)
  )
cat("  Samples:", nrow(sample_annot), "\n")
cat("  Age groups:", paste(names(table(sample_annot$AgeRange)), collapse = ", "), "\n")
cat("  Sample types:", paste(names(table(sample_annot$Group)), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Step 4: Prepare DESeq2 Input
# -----------------------------------------------------------------------------
cat("Step 4: Preparing DESeq2 input...\n")

# Transpose and join
count_data_t <- count_data_filtered %>%
  tibble::column_to_rownames("GeneSymb") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleName")

count_annot <- dplyr::inner_join(
  sample_annot[, c("SampleName", "SampleNameGroup", "AgeRange", "Group")],
  count_data_t,
  by = "SampleName"
) %>%
  dplyr::select(-SampleName)

# Create count matrix (genes x samples)
count_matrix <- count_annot %>%
  tibble::column_to_rownames("SampleNameGroup") %>%
  dplyr::select(-AgeRange, -Group) %>%
  t() %>%
  as.data.frame()

# Create metadata
metadata <- data.frame(
  SampleNameGroup = colnames(count_matrix),
  row.names = colnames(count_matrix)
)
metadata <- metadata %>%
  dplyr::mutate(
    AgeRange = gsub(".*_(.*)Tumor.*", "\\1", SampleNameGroup),
    Group = ifelse(grepl("TumorAdj", SampleNameGroup), "TumorAdj", "Tumor")
  )

cat("  Count matrix:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# -----------------------------------------------------------------------------
# Step 5: Create DESeq2 Object with CORRECTED Design
# -----------------------------------------------------------------------------
cat("Step 5: Creating DESeq2 object...\n")

# BIOSTATISTICAL FIX: Use proper design formula instead of ~1
# This accounts for AgeRange and Group effects during normalization
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_matrix),
  colData = metadata,
  design = ~ AgeRange + Group  # CORRECTED: was design = ~1
)

cat("  DESeq2 object created with design: ~ AgeRange + Group\n")

# -----------------------------------------------------------------------------
# Step 6: Variance Stabilizing Transformation
# -----------------------------------------------------------------------------
cat("Step 6: Applying variance stabilizing transformation...\n")
vst_data <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_matrix <- assay(vst_data) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneSymb")

# Remove duplicate genes
vst_matrix_nodup <- vst_matrix %>%
  dplyr::filter(!duplicated(GeneSymb))
cat("  VST matrix:", nrow(vst_matrix_nodup), "genes\n")

# -----------------------------------------------------------------------------
# Step 7: Generate QC Plots
# -----------------------------------------------------------------------------
cat("Step 7: Generating QC plots...\n")

pdf(file.path(output_dir, "qc_plots.pdf"), width = 12, height = 10)

# PCA plot
pca_data <- plotPCA(vst_data, intgroup = c("AgeRange", "Group"), returnData = TRUE)
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = AgeRange, shape = Group)) +
  geom_point(size = 3) +
  theme_bw() +
  ggtitle("PCA: Colored by Age Range, Shaped by Sample Type")
print(pca_plot)

# Sample distance heatmap
sample_dists <- dist(t(assay(vst_data)))
sample_dist_matrix <- as.matrix(sample_dists)
pheatmap::pheatmap(
  sample_dist_matrix,
  main = "Sample Distance Matrix",
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists
)

# Library size distribution
lib_sizes <- colSums(counts(dds))
lib_df <- data.frame(
  Sample = names(lib_sizes),
  LibrarySize = lib_sizes,
  Group = metadata$Group,
  AgeRange = metadata$AgeRange
)
lib_plot <- ggplot(lib_df, aes(x = reorder(Sample, LibrarySize), y = LibrarySize / 1e6, fill = Group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() +
  labs(x = "Sample", y = "Library Size (millions)", title = "Library Size Distribution") +
  theme(axis.text.y = element_text(size = 6))
print(lib_plot)

dev.off()
cat("  QC plots saved to:", file.path(output_dir, "qc_plots.pdf"), "\n")

# -----------------------------------------------------------------------------
# Step 8: Save Outputs
# -----------------------------------------------------------------------------
cat("Step 8: Saving outputs...\n")
saveRDS(dds, file.path(output_dir, "dds_norm_AllAgeGroup.rds"))
saveRDS(vst_matrix_nodup, file.path(output_dir, "vst_normalized_matrix.rds"))
saveRDS(sample_annot, file.path(output_dir, "sample_annotation.rds"))

# Memory cleanup - remove large intermediate objects
rm(count_data, count_data_t, count_annot, gene_annot, vst_data)
gc()

cat("\n=== Preprocessing complete ===\n")
cat("Outputs saved to:", output_dir, "\n")
