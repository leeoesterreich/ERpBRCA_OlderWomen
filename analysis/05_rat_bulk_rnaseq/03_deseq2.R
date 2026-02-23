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
#   - analysis/05_rat_bulk_rnaseq/outputs/normalized_tpm.csv (for PAM50)

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
cat("Count dir:", count_dir, "\n")
cat("Data dir:", data_dir, "\n")

# -----------------------------------------------------------------------------
# Step 1: Load count data
# -----------------------------------------------------------------------------
cat("\nStep 1: Loading count data...\n")

count_files <- list.files(count_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(count_files) == 0) {
  stop("No count files found in: ", count_dir)
}

count_list <- lapply(count_files, function(f) {
  df <- fread(f, header = FALSE, col.names = c("gene_id", "count"))
  df$sample <- gsub("\\.txt$", "", basename(f))
  df
})

counts_long <- rbindlist(count_list)
counts_wide <- dcast(counts_long, gene_id ~ sample, value.var = "count")

# Remove HTSeq summary lines
counts_wide <- counts_wide[!grepl("^__", gene_id)]
countdata <- as.matrix(counts_wide[, -1, with = FALSE])
rownames(countdata) <- counts_wide$gene_id

cat("  Count matrix:", nrow(countdata), "genes x", ncol(countdata), "samples\n")

# -----------------------------------------------------------------------------
# Step 2: Load sample metadata
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading sample metadata...\n")

metadata_file <- file.path(data_dir, "sample_metadata.csv")
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

coldata <- read.csv(metadata_file)
rownames(coldata) <- coldata$SAMPLE
coldata$TYPE <- factor(coldata$TYPE, levels = c("CONTROL", "TEST"))

# Match sample order
common_samples <- intersect(colnames(countdata), rownames(coldata))
if (length(common_samples) == 0) {
  cat("  Count file samples:", paste(colnames(countdata), collapse = ", "), "\n")
  cat("  Metadata samples:", paste(rownames(coldata), collapse = ", "), "\n")
  stop("No matching samples between counts and metadata!")
}

countdata <- countdata[, common_samples]
coldata <- coldata[common_samples, ]

cat("  Samples:", nrow(coldata), "\n")
cat("  Groups:", paste(table(coldata$TYPE), collapse = " vs "), "\n")

# -----------------------------------------------------------------------------
# Step 3: Run DESeq2
# -----------------------------------------------------------------------------
cat("\nStep 3: Running DESeq2...\n")

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ TYPE
)

dds <- DESeq(dds)

# Get results (TEST vs CONTROL)
results <- results(dds, contrast = c("TYPE", "TEST", "CONTROL"))

cat("  Total genes tested:", sum(!is.na(results$padj)), "\n")

# -----------------------------------------------------------------------------
# Step 4: Add gene symbols
# -----------------------------------------------------------------------------
cat("\nStep 4: Adding gene symbols...\n")

res_df <- as.data.frame(results) %>%
  tibble::rownames_to_column("ensembl_gene_id")

# Add normalized counts
norm_counts <- counts(dds, normalized = TRUE)
res_df <- merge(res_df,
                as.data.frame(norm_counts) %>% tibble::rownames_to_column("ensembl_gene_id"),
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
  cat("  Gene symbols added\n")
}, error = function(e) {
  cat("  Warning: Could not connect to biomaRt:", conditionMessage(e), "\n")
  res_df$external_gene_name <<- NA
})

# -----------------------------------------------------------------------------
# Step 5: Save results
# -----------------------------------------------------------------------------
cat("\nStep 5: Saving results...\n")

# Save all results
write.csv(res_df, file.path(output_dir, "deseq2_results.csv"), row.names = FALSE)
cat("  Saved: deseq2_results.csv\n")

# Filter by FDR < 0.05
res_sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

write.csv(res_sig, file.path(output_dir, "deseq2_results_significant.csv"), row.names = FALSE)
cat("  Significant genes (FDR < 0.05):", nrow(res_sig), "\n")

# Summary
cat("\nSummary:\n")
cat("  Upregulated (log2FC > 0):", sum(res_sig$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Downregulated (log2FC < 0):", sum(res_sig$log2FoldChange < 0, na.rm = TRUE), "\n")

# Save normalized counts for PAM50
norm_tpm <- sweep(norm_counts, 2, colSums(norm_counts), "/") * 1e6
write.csv(norm_tpm, file.path(output_dir, "normalized_tpm.csv"))
cat("  Saved: normalized_tpm.csv\n")

cat("\n=== DESeq2 complete ===\n")
