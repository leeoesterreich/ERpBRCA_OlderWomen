#!/usr/bin/env Rscript
# analysis/01_human_bulk_rnaseq/01_preprocess.R
# Preprocess human bulk RNA-seq data: read counts/TPM, filter genes, normalize
#
# Inputs (Count mode - preferred):
#   - data/human_bulk_rnaseq/raw/HumanERpAge_39404g168s_FeatureCount.txt
#   - data/human_bulk_rnaseq/raw/HumanERpAge_BulkRNAseq_SampleInformation.txt
#   - data/human_bulk_rnaseq/external/Homo_sapiens.gene_info.txt
#
# Inputs (TPM mode - fallback when count data unavailable):
#   - data/human_bulk_rnaseq/raw/HumanERpAge_39404g168s_TPMlog2.txt
#   - data/human_bulk_rnaseq/raw/HumanERpAge_BulkRNAseq_SampleInformation.txt
#
# Outputs:
#   - analysis/01_human_bulk_rnaseq/outputs/dds_norm_AllAgeGroup.rds (count mode only)
#   - analysis/01_human_bulk_rnaseq/outputs/vst_normalized_matrix.rds
#   - analysis/01_human_bulk_rnaseq/outputs/qc_plots.pdf

# Set random seed for reproducibility
set.seed(12345)

# Load libraries
suppressPackageStartupMessages({
  library(tidyr)
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

# Define paths - use commandArgs to get script directory when run via Rscript
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "../.."))
data_dir <- file.path(project_root, "data/human_bulk_rnaseq")
output_dir <- file.path(script_dir, "outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input files
count_file <- file.path(data_dir, "raw/HumanERpAge_39404g168s_FeatureCount.txt")
tpm_file <- file.path(data_dir, "raw/HumanERpAge_39404g168s_TPMlog2.txt")
sample_file <- file.path(data_dir, "raw/HumanERpAge_BulkRNAseq_SampleInformation.txt")
gene_annot_file <- file.path(data_dir, "external/Homo_sapiens.gene_info.txt")

cat("=== Human Bulk RNA-seq Preprocessing ===\n")
cat("Project root:", project_root, "\n")
cat("Output directory:", output_dir, "\n\n")

# Verify input files exist - use TPM if count data not available
cat("Checking input files...\n")
use_tpm_only <- FALSE
if (!file.exists(count_file)) {
  cat("  Note: Count data not found, will use TPM data directly\n")
  use_tpm_only <- TRUE
  check_file_exists(tpm_file, "TPM data file")
} else {
  check_file_exists(count_file, "Count data file")
}
check_file_exists(sample_file, "Sample annotation file")

# Gene annotation is optional when using TPM
has_gene_annot <- file.exists(gene_annot_file)
if (has_gene_annot) {
  check_file_exists(gene_annot_file, "Gene annotation file")
} else {
  cat("  Note: Gene annotation file not found, skipping protein-coding filter\n")
}
cat("\n")

# =============================================================================
# TPM-ONLY MODE: Simplified preprocessing when only TPM data available
# =============================================================================
if (use_tpm_only) {
  cat("*** Running in TPM-only mode (no DESeq2 normalization) ***\n\n")

  # -------------------------------------------------------------------------
  # Step 1: Read Sample Annotation
  # -------------------------------------------------------------------------
  cat("Step 1: Reading sample annotation...\n")
  sample_annot <- fread(sample_file, stringsAsFactors = FALSE, header = TRUE)

  # Handle column name variations
  if ("SampleID" %in% colnames(sample_annot)) {
    colnames(sample_annot)[colnames(sample_annot) == "SampleID"] <- "SampleName"
  } else if (ncol(sample_annot) >= 2 && colnames(sample_annot)[1] != "SampleName") {
    colnames(sample_annot)[1] <- "SampleName"
  }

  # Map TissueType to Group for consistency with downstream scripts
  if ("TissueType" %in% colnames(sample_annot) && !"Group" %in% colnames(sample_annot)) {
    sample_annot$Group <- gsub("Tumor_Adj", "TumorAdj", sample_annot$TissueType)
    sample_annot$Group <- gsub("^Tumor$", "Tumor", sample_annot$Group)
  }

  # Create Age alias from ChronologicalAge for downstream scripts
  if ("ChronologicalAge" %in% colnames(sample_annot) && !"Age" %in% colnames(sample_annot)) {
    sample_annot$Age <- sample_annot$ChronologicalAge
  }

  sample_annot <- sample_annot %>%
    dplyr::mutate(
      AgeRange_Group = paste0(AgeRange, "_", Group),
      SampleNameGroup = paste0(SampleName, "_", AgeRange, Group)
    )
  cat("  Samples:", nrow(sample_annot), "\n")
  cat("  Age groups:", paste(names(table(sample_annot$AgeRange)), collapse = ", "), "\n")
  cat("  Sample types:", paste(names(table(sample_annot$Group)), collapse = ", "), "\n")

  # -------------------------------------------------------------------------
  # Step 2: Read TPM Data
  # -------------------------------------------------------------------------
  cat("Step 2: Reading TPM data (log2-transformed)...\n")
  tpm_data <- fread(tpm_file, header = TRUE, stringsAsFactors = FALSE)
  colnames(tpm_data)[1] <- "GeneSymb"

  # Strip "_LEE*" suffix from sample names to match annotation
  colnames(tpm_data) <- gsub("_LEE.*", "", colnames(tpm_data))

  cat("  Raw genes:", nrow(tpm_data), "\n")
  cat("  Samples in data:", ncol(tpm_data) - 1, "\n")

  # -------------------------------------------------------------------------
  # Step 3: Filter to Protein-Coding Genes (optional)
  # -------------------------------------------------------------------------
  if (has_gene_annot) {
    cat("Step 3: Filtering to protein-coding genes...\n")
    gene_annot <- fread(gene_annot_file, header = TRUE, stringsAsFactors = FALSE)
    gene_annot_filtered <- gene_annot %>%
      dplyr::select(Symbol, type_of_gene) %>%
      dplyr::filter(type_of_gene == "protein-coding", !grepl("^LOC\\d+", Symbol))

    n_before <- nrow(tpm_data)
    tpm_data <- tpm_data %>%
      dplyr::filter(GeneSymb %in% gene_annot_filtered$Symbol)
    cat("  Protein-coding genes retained:", nrow(tpm_data), "\n")
    cat("  Genes removed:", n_before - nrow(tpm_data), "\n")
  } else {
    cat("Step 3: Skipping gene filter (no annotation)\n")
  }

  # -------------------------------------------------------------------------
  # Step 4: Match Samples and Create Output Matrix
  # -------------------------------------------------------------------------
  cat("Step 4: Matching samples to annotation...\n")

  # Get sample names from TPM data (excluding GeneSymb column)
  tpm_sample_names <- colnames(tpm_data)[-1]

  # Match with annotation
  matched_samples <- sample_annot %>%
    dplyr::filter(SampleName %in% tpm_sample_names)

  cat("  Samples matched:", nrow(matched_samples), "of", nrow(sample_annot), "\n")

  if (nrow(matched_samples) == 0) {
    stop("ERROR: No samples matched between TPM data and annotation")
  }

  # Create normalized matrix with proper column names
  tpm_matrix <- tpm_data %>%
    dplyr::select(GeneSymb, all_of(matched_samples$SampleName))

  # Rename columns to SampleNameGroup format
  new_colnames <- c("GeneSymb", matched_samples$SampleNameGroup)
  colnames(tpm_matrix) <- new_colnames

  # Remove duplicate genes
  tpm_matrix_nodup <- tpm_matrix %>%
    dplyr::filter(!duplicated(GeneSymb))

  cat("  Final matrix:", nrow(tpm_matrix_nodup), "genes x", ncol(tpm_matrix_nodup) - 1, "samples\n")

  # -------------------------------------------------------------------------
  # Step 5: Create Metadata
  # -------------------------------------------------------------------------
  cat("Step 5: Creating metadata...\n")
  metadata <- data.frame(
    SampleNameGroup = matched_samples$SampleNameGroup,
    SampleName = matched_samples$SampleName,
    AgeRange = matched_samples$AgeRange,
    Group = matched_samples$Group,
    row.names = matched_samples$SampleNameGroup
  )

  # -------------------------------------------------------------------------
  # Step 6: Generate QC Plots
  # -------------------------------------------------------------------------
  cat("Step 6: Generating QC plots...\n")

  # Prepare numeric matrix for PCA
  expr_matrix <- tpm_matrix_nodup %>%
    tibble::column_to_rownames("GeneSymb") %>%
    as.matrix()

  # Run PCA
  pca_result <- prcomp(t(expr_matrix), center = TRUE, scale. = TRUE)
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    SampleNameGroup = rownames(pca_result$x)
  )
  pca_df <- dplyr::left_join(pca_df, metadata, by = "SampleNameGroup")

  # Calculate variance explained
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

  pdf(file.path(output_dir, "qc_plots.pdf"), width = 12, height = 10)

  # PCA plot
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = AgeRange, shape = Group)) +
    geom_point(size = 3) +
    theme_bw() +
    labs(
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      title = "PCA: Colored by Age Range, Shaped by Sample Type",
      subtitle = "Based on log2(TPM) expression"
    )
  print(pca_plot)

  # Sample correlation heatmap
  sample_cor <- cor(expr_matrix)
  pheatmap::pheatmap(
    sample_cor,
    main = "Sample Correlation Matrix (log2 TPM)",
    annotation_col = metadata[, c("AgeRange", "Group"), drop = FALSE]
  )

  # Expression distribution per sample
  expr_long <- tpm_matrix_nodup %>%
    tidyr::pivot_longer(-GeneSymb, names_to = "Sample", values_to = "log2TPM") %>%
    dplyr::left_join(metadata, by = c("Sample" = "SampleNameGroup"))

  dist_plot <- ggplot(expr_long, aes(x = log2TPM, fill = AgeRange)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Group) +
    theme_bw() +
    labs(title = "Expression Distribution by Age and Tissue Type", x = "log2(TPM)")
  print(dist_plot)

  dev.off()
  cat("  QC plots saved to:", file.path(output_dir, "qc_plots.pdf"), "\n")

  # -------------------------------------------------------------------------
  # Step 7: Save Outputs
  # -------------------------------------------------------------------------
  cat("Step 7: Saving outputs...\n")

  # Save the TPM matrix as the normalized matrix (already log2 transformed)
  saveRDS(tpm_matrix_nodup, file.path(output_dir, "vst_normalized_matrix.rds"))
  saveRDS(matched_samples, file.path(output_dir, "sample_annotation.rds"))
  saveRDS(metadata, file.path(output_dir, "metadata.rds"))

  # Note: No dds object when using TPM-only mode
  cat("  Note: No DESeq2 object saved (TPM-only mode)\n")

  cat("\n=== Preprocessing complete (TPM-only mode) ===\n")
  cat("Outputs saved to:", output_dir, "\n")

} else {
  # ===========================================================================
  # COUNT MODE: Full DESeq2 preprocessing
  # ===========================================================================

  # Load DESeq2 only when needed
  suppressPackageStartupMessages(library(DESeq2))

  # -------------------------------------------------------------------------
  # Step 1: Read Gene Annotation
  # -------------------------------------------------------------------------
  cat("Step 1: Reading gene annotation...\n")
  gene_annot <- fread(gene_annot_file, header = TRUE, stringsAsFactors = FALSE)
  n_genes_raw <- nrow(gene_annot)
  gene_annot_filtered <- gene_annot %>%
    dplyr::select(Symbol, type_of_gene) %>%
    dplyr::filter(type_of_gene == "protein-coding", !grepl("^LOC\\d+", Symbol))
  cat("  Total genes in annotation:", n_genes_raw, "\n")
  cat("  Protein-coding genes (after filter):", nrow(gene_annot_filtered), "\n")
  cat("  Genes removed:", n_genes_raw - nrow(gene_annot_filtered), "\n")

  # -------------------------------------------------------------------------
  # Step 2: Read Count Data
  # -------------------------------------------------------------------------
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

  # -------------------------------------------------------------------------
  # Step 3: Read Sample Annotation
  # -------------------------------------------------------------------------
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

  # -------------------------------------------------------------------------
  # Step 4: Prepare DESeq2 Input
  # -------------------------------------------------------------------------
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

  # -------------------------------------------------------------------------
  # Step 5: Create DESeq2 Object with CORRECTED Design
  # -------------------------------------------------------------------------
  cat("Step 5: Creating DESeq2 object...\n")

  # BIOSTATISTICAL FIX: Use proper design formula instead of ~1
  # This accounts for AgeRange and Group effects during normalization
  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(count_matrix),
    colData = metadata,
    design = ~ AgeRange + Group  # CORRECTED: was design = ~1
  )

  cat("  DESeq2 object created with design: ~ AgeRange + Group\n")

  # -------------------------------------------------------------------------
  # Step 6: Variance Stabilizing Transformation
  # -------------------------------------------------------------------------
  cat("Step 6: Applying variance stabilizing transformation...\n")
  vst_data <- varianceStabilizingTransformation(dds, blind = FALSE)
  vst_matrix <- assay(vst_data) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("GeneSymb")

  # Remove duplicate genes
  vst_matrix_nodup <- vst_matrix %>%
    dplyr::filter(!duplicated(GeneSymb))
  cat("  VST matrix:", nrow(vst_matrix_nodup), "genes\n")

  # -------------------------------------------------------------------------
  # Step 7: Generate QC Plots
  # -------------------------------------------------------------------------
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

  # -------------------------------------------------------------------------
  # Step 8: Save Outputs
  # -------------------------------------------------------------------------
  cat("Step 8: Saving outputs...\n")
  saveRDS(dds, file.path(output_dir, "dds_norm_AllAgeGroup.rds"))
  saveRDS(vst_matrix_nodup, file.path(output_dir, "vst_normalized_matrix.rds"))
  saveRDS(sample_annot, file.path(output_dir, "sample_annotation.rds"))

  # Memory cleanup - remove large intermediate objects
  rm(count_data, count_data_t, count_annot, gene_annot, vst_data)
  gc()

  cat("\n=== Preprocessing complete ===\n")
  cat("Outputs saved to:", output_dir, "\n")
}
