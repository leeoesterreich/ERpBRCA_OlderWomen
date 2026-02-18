#!/usr/bin/env Rscript
# analysis/05_rat_bulk_rnaseq/03_deseq2.R
# Differential expression analysis for rat bulk RNA-seq
#
# Inputs:
#   - analysis/05_rat_bulk_rnaseq/outputs/counts/*.txt (HTSeq counts)
#   - data/rat_bulk_rnaseq/sample_metadata.csv
#
# Outputs:
#   - analysis/05_rat_bulk_rnaseq/outputs/deseq2_results.csv
#   - analysis/05_rat_bulk_rnaseq/outputs/deseq2_results_significant.csv

set.seed(12345)

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(biomaRt)
  library(dplyr)
})

# Define paths
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
output_dir <- file.path(script_dir, "outputs")
count_dir <- file.path(output_dir, "counts")
data_dir <- file.path(project_root, "data/rat_bulk_rnaseq")

cat("=== DESeq2 Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load count data
# -----------------------------------------------------------------------------
cat("Step 1: Loading count data...\n")

# Read and combine count files
count_files <- list.files(count_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(count_files) == 0) {
  stop("No count files found in: ", count_dir)
}

count_list <- lapply(count_files, function(f) {
  df <- fread(f, header = FALSE, col.names = c("gene_id", basename(f)))
  df
})

# Merge all count files
countdata <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), count_list)
countdata <- countdata %>%
  filter(!grepl("^__", gene_id)) %>%  # Remove HTSeq summary lines
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

# Clean column names
colnames(countdata) <- gsub("\\.txt$", "", colnames(countdata))

cat("  Count matrix:", nrow(countdata), "genes x", ncol(countdata), "samples\n")

# -----------------------------------------------------------------------------
# Step 2: Load sample metadata
# -----------------------------------------------------------------------------
cat("Step 2: Loading sample metadata...\n")

metadata_file <- file.path(data_dir, "sample_metadata.csv")
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

coldata <- read.csv(metadata_file, row.names = 1)
coldata$TYPE <- factor(coldata$TYPE)

# Ensure sample order matches
countdata <- countdata[, rownames(coldata)]

cat("  Samples:", nrow(coldata), "\n")
cat("  Groups:", paste(levels(coldata$TYPE), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# Step 3: Run DESeq2
# -----------------------------------------------------------------------------
cat("Step 3: Running DESeq2...\n")

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ TYPE
)

dds <- DESeq(dds)

# Get results with contrast
results <- results(dds, contrast = c("TYPE", "TEST", "CONTROL"))

cat("  Total genes tested:", nrow(results), "\n")

# -----------------------------------------------------------------------------
# Step 4: Add gene symbols and apply FDR correction
# -----------------------------------------------------------------------------
cat("Step 4: Adding gene symbols...\n")

res_df <- as.data.frame(results) %>%
  tibble::rownames_to_column("ensembl_gene_id")

# Add normalized counts
norm_counts <- counts(dds, normalized = TRUE)
res_df <- merge(res_df, as.data.frame(norm_counts) %>% tibble::rownames_to_column("ensembl_gene_id"),
                by = "ensembl_gene_id", all.x = TRUE)

# Get gene symbols from biomaRt
tryCatch({
  mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  gene_symbols <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "external_gene_name"),
    values = res_df$ensembl_gene_id,
    mart = mart
  )
  res_df <- merge(res_df, gene_symbols, by = "ensembl_gene_id", all.x = TRUE)
}, error = function(e) {
  cat("  Warning: Could not connect to biomaRt, skipping gene symbol annotation\n")
  res_df$external_gene_name <- NA
})

# -----------------------------------------------------------------------------
# Step 5: Save results with FDR filtering
# -----------------------------------------------------------------------------
cat("Step 5: Saving results...\n")

# Save all results
write.csv(res_df, file.path(output_dir, "deseq2_results.csv"), row.names = FALSE)
cat("  Saved: deseq2_results.csv\n")

# FIX: Filter by adjusted p-value (FDR < 0.05), not raw p-value
res_sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

write.csv(res_sig, file.path(output_dir, "deseq2_results_significant.csv"), row.names = FALSE)
cat("  Significant genes (FDR < 0.05):", nrow(res_sig), "\n")
cat("  Saved: deseq2_results_significant.csv\n")

# Summary statistics
cat("\nSummary:\n")
cat("  Upregulated (FDR < 0.05, log2FC > 0):", sum(res_sig$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Downregulated (FDR < 0.05, log2FC < 0):", sum(res_sig$log2FoldChange < 0, na.rm = TRUE), "\n")

cat("\n=== DESeq2 complete ===\n")
