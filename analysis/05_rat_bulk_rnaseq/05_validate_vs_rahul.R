#!/usr/bin/env Rscript
# analysis/05_rat_bulk_rnaseq/05_validate_vs_rahul.R
# Compare pipeline outputs to Rahul's original analysis

set.seed(12345)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

cat("=== Validation Against Rahul's Results ===\n\n")

# Paths
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
rahul_dir <- "/ix1/alee/LO_LAB/Personal/Rahul/Neil_RNAseq"

# Sample mapping (our names to Rahul's)
sample_map <- c(
  "102-FF-Tumor" = "102-FF",
  "107-L4t2-Tumor" = "107-L4t2",
  "116-R4-Tumor" = "116-R4-Tumor",
  "157-R2-Tumor" = "157-R2-Tumor",
  "158-L2-Tumor" = "158-L2-Tumor",
  "167-L2-Tumor" = "167-L2-Tumor"
)

results <- list()

# -----------------------------------------------------------------------------
# Check 1: HTSeq counts correlation
# -----------------------------------------------------------------------------
cat("Check 1: HTSeq Counts Correlation\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

our_count_dir <- file.path(output_dir, "counts")
rahul_count_dir <- file.path(rahul_dir, "5_Count_file")

correlations <- sapply(names(sample_map), function(our_name) {
  rahul_name <- sample_map[our_name]

  our_file <- file.path(our_count_dir, paste0(our_name, ".txt"))
  rahul_file <- file.path(rahul_count_dir, paste0(rahul_name, ".txt"))

  if (!file.exists(our_file)) {
    cat(sprintf("  %s: OUR FILE MISSING\n", our_name))
    return(NA)
  }
  if (!file.exists(rahul_file)) {
    cat(sprintf("  %s: RAHUL FILE MISSING\n", our_name))
    return(NA)
  }

  our_counts <- fread(our_file, header = FALSE, col.names = c("gene", "count"))
  rahul_counts <- fread(rahul_file, header = FALSE, col.names = c("gene", "count"))

  # Merge on gene ID
  merged <- merge(our_counts, rahul_counts, by = "gene", suffixes = c("_ours", "_rahul"))
  merged <- merged[!grepl("^__", gene)]  # Remove summary rows

  # Calculate correlation
  r <- cor(merged$count_ours, merged$count_rahul, method = "pearson")
  cat(sprintf("  %s: r = %.6f %s\n", our_name, r, ifelse(r >= 0.99, "PASS", "FAIL")))
  return(r)
})

results$count_correlations <- correlations
results$count_pass <- all(correlations >= 0.99, na.rm = TRUE)
cat(sprintf("\nOverall: %s (min r = %.4f)\n\n",
            ifelse(results$count_pass, "PASS", "FAIL"),
            min(correlations, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# Check 2: DESeq2 significant gene overlap
# -----------------------------------------------------------------------------
cat("Check 2: DESeq2 Significant Gene Overlap\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

our_de <- tryCatch({
  fread(file.path(output_dir, "deseq2_results_significant.csv"))
}, error = function(e) NULL)

rahul_de <- tryCatch({
  fread(file.path(rahul_dir, "6_DEseq2/DESeq2_results_with_symbols.csv"))
}, error = function(e) NULL)

if (!is.null(our_de) && !is.null(rahul_de)) {
  # Filter Rahul's to FDR < 0.05
  rahul_sig <- rahul_de[!is.na(padj) & padj < 0.05]

  our_genes <- our_de$ensembl_gene_id
  rahul_genes <- rahul_sig$Row.names

  overlap <- length(intersect(our_genes, rahul_genes))
  union_size <- length(union(our_genes, rahul_genes))
  jaccard <- overlap / union_size

  cat(sprintf("  Our significant genes: %d\n", length(our_genes)))
  cat(sprintf("  Rahul's significant genes: %d\n", length(rahul_genes)))
  cat(sprintf("  Overlap: %d (%.1f%%)\n", overlap, 100 * overlap / min(length(our_genes), length(rahul_genes))))
  cat(sprintf("  Jaccard index: %.3f\n", jaccard))

  results$de_overlap <- overlap
  results$de_jaccard <- jaccard
  results$de_pass <- (overlap / min(length(our_genes), length(rahul_genes))) >= 0.90
  cat(sprintf("\nOverall: %s\n\n", ifelse(results$de_pass, "PASS", "FAIL")))
} else {
  cat("  Could not load DESeq2 results\n\n")
  results$de_pass <- NA
}

# -----------------------------------------------------------------------------
# Check 3: PAM50 subtype concordance
# -----------------------------------------------------------------------------
cat("Check 3: PAM50 Subtype Concordance\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

our_pam50 <- tryCatch({
  fread(file.path(output_dir, "pam50_subtypes.csv"))
}, error = function(e) NULL)

rahul_pam50 <- tryCatch({
  fread(file.path(rahul_dir, "PAM50/PAM50.csv"))
}, error = function(e) NULL)

if (!is.null(our_pam50) && !is.null(rahul_pam50)) {
  # Compare subtypes
  cat("  Sample comparisons:\n")
  matches <- 0
  total <- 0

  for (i in 1:nrow(our_pam50)) {
    our_sample <- our_pam50$Sample[i]
    our_subtype <- our_pam50$Subtype[i]

    # Find matching sample in Rahul's results
    rahul_row <- rahul_pam50[grepl(gsub("-Tumor", "", our_sample), rahul_pam50[[1]], ignore.case = TRUE)]

    if (nrow(rahul_row) > 0) {
      rahul_subtype <- rahul_row$Subtype[1]
      match <- our_subtype == rahul_subtype
      matches <- matches + match
      total <- total + 1
      cat(sprintf("    %s: %s vs %s %s\n",
                  our_sample, our_subtype, rahul_subtype,
                  ifelse(match, "MATCH", "MISMATCH")))
    }
  }

  results$pam50_matches <- matches
  results$pam50_total <- total
  results$pam50_pass <- matches == total
  cat(sprintf("\nOverall: %s (%d/%d match)\n\n",
              ifelse(results$pam50_pass, "PASS", "FAIL"), matches, total))
} else {
  cat("  Could not load PAM50 results\n\n")
  results$pam50_pass <- NA
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("=== VALIDATION SUMMARY ===\n")
cat(sprintf("HTSeq Counts:  %s\n", ifelse(isTRUE(results$count_pass), "PASS", "FAIL")))
cat(sprintf("DESeq2 Genes:  %s\n", ifelse(isTRUE(results$de_pass), "PASS", "FAIL")))
cat(sprintf("PAM50 Types:   %s\n", ifelse(isTRUE(results$pam50_pass), "PASS", "FAIL")))

all_pass <- isTRUE(results$count_pass) && isTRUE(results$de_pass) && isTRUE(results$pam50_pass)
cat(sprintf("\nOVERALL: %s\n", ifelse(all_pass, "ALL CHECKS PASSED", "SOME CHECKS FAILED")))

# Save results
saveRDS(results, file.path(output_dir, "validation_results.rds"))
cat("\nResults saved to: validation_results.rds\n")
