#!/usr/bin/env Rscript
# analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R
# PAM50 molecular subtyping for rat tumors
#
# Inputs:
#   - analysis/05_rat_bulk_rnaseq/outputs/counts/*.txt (HTSeq counts)
#   - Human-to-rat PAM50 gene ortholog mapping
#
# Outputs:
#   - analysis/05_rat_bulk_rnaseq/outputs/pam50_subtypes.csv
#   - analysis/05_rat_bulk_rnaseq/outputs/pam50_heatmap.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(genefu)
  library(dplyr)
  library(pheatmap)
  library(data.table)
  library(biomaRt)
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
output_dir <- file.path(script_dir, "outputs")

cat("=== PAM50 Subtyping ===\n")

# Load PAM50 reference data
data("pam50")
data("pam50.robust")

# -----------------------------------------------------------------------------
# Step 1: Define human-to-rat ortholog mapping
# -----------------------------------------------------------------------------
cat("Step 1: Setting up gene ortholog mapping...\n")

pam50_genes_human <- rownames(pam50$centroids)

human_to_rat <- c(
  "ACTR3B" = "Actr3b", "ANLN" = "Anln", "BAG1" = "Bag1", "BCL2" = "Bcl2",
  "BIRC5" = "Birc5", "BLVRA" = "Blvra", "CCNB1" = "Ccnb1", "CCNE1" = "Ccne1",
  "CDC20" = "Cdc20", "CDCA1" = "Nuf2", "NUF2" = "Nuf2", "CDC6" = "Cdc6",
  "CDH3" = "Cdh3", "CENPF" = "Cenpf", "CEP55" = "Cep55", "CXXC5" = "Cxxc5",
  "EGFR" = "Egfr", "ERBB2" = "Erbb2", "ESR1" = "Esr1", "EXO1" = "Exo1",
  "FGFR4" = "Fgfr4", "FOXA1" = "Foxa1", "FOXC1" = "Foxc1", "GPR160" = "Gpr160",
  "GRB7" = "Grb7", "KIF2C" = "Kif2c", "KRT14" = "Krt14", "KRT17" = "Krt17",
  "KRT5" = "Krt5", "MAPT" = "Mapt", "MDM2" = "Mdm2", "MELK" = "Melk",
  "MIA" = "Mia", "MKI67" = "Mki67", "MLPH" = "Mlph", "MMP11" = "Mmp11",
  "MYBL2" = "Mybl2", "MYC" = "Myc", "NAT1" = "Nat1", "ORC6" = "Orc6",
  "PGR" = "Pgr", "PHGDH" = "Phgdh", "PTTG1" = "Pttg1", "RRM2" = "Rrm2",
  "SFRP1" = "Sfrp1", "SLC39A6" = "Slc39a6", "TMEM45B" = "Tmem45b", "TYMS" = "Tyms",
  "UBE2C" = "Ube2c", "UBE2T" = "Ube2t",
  # Alternative gene names
  "KNTC2" = "Ndc80",
  "ORC6L" = "Orc6"
)

pam50_genes_rat <- human_to_rat[pam50_genes_human]

cat("  PAM50 genes mapped:", sum(!is.na(pam50_genes_rat)), "/", length(pam50_genes_human), "\n")

# -----------------------------------------------------------------------------
# Step 2: Load expression data
# -----------------------------------------------------------------------------
cat("Step 2: Loading expression data...\n")

# Assume normalized TPM data exists or compute from counts
tpm_file <- file.path(output_dir, "normalized_tpm.csv")
if (file.exists(tpm_file)) {
  pam50_data <- read.csv(tpm_file, row.names = 1)
} else {
  stop("Normalized TPM file not found. Run normalization first: ", tpm_file)
}

cat("  Expression matrix:", nrow(pam50_data), "genes x", ncol(pam50_data), "samples\n")

# -----------------------------------------------------------------------------
# Step 2b: Convert Ensembl IDs to gene symbols
# -----------------------------------------------------------------------------
cat("Step 2b: Converting Ensembl IDs to gene symbols...\n")

# Check if row names are Ensembl IDs
if (grepl("^ENSRNOG", rownames(pam50_data)[1])) {
  cat("  Detected Ensembl IDs, fetching gene symbols from biomaRt...\n")

  mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  gene_map <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "external_gene_name"),
    values = rownames(pam50_data),
    mart = mart
  )

  # Create mapping, keep only unique mappings
  gene_map <- gene_map[gene_map$external_gene_name != "", ]
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]

  # Map Ensembl IDs to gene symbols
  id_to_symbol <- setNames(gene_map$external_gene_name, gene_map$ensembl_gene_id)
  new_rownames <- id_to_symbol[rownames(pam50_data)]

  # Keep rows that have valid gene symbols
  valid_idx <- !is.na(new_rownames) & new_rownames != ""
  pam50_data <- pam50_data[valid_idx, ]
  new_rownames <- new_rownames[valid_idx]

  # Handle duplicates BEFORE setting row names: keep highest expressing version
  # For duplicated gene symbols, keep the one with highest mean expression
  dup_symbols <- unique(new_rownames[duplicated(new_rownames)])
  keep_idx <- rep(TRUE, nrow(pam50_data))

  for (sym in dup_symbols) {
    dup_idx <- which(new_rownames == sym)
    mean_expr <- rowMeans(pam50_data[dup_idx, , drop = FALSE])
    keep_one <- dup_idx[which.max(mean_expr)]
    keep_idx[dup_idx] <- FALSE
    keep_idx[keep_one] <- TRUE
  }

  pam50_data <- pam50_data[keep_idx, ]
  new_rownames <- new_rownames[keep_idx]
  rownames(pam50_data) <- new_rownames

  cat("  Converted:", nrow(pam50_data), "genes with unique symbols\n")
}

# -----------------------------------------------------------------------------
# Step 3: Check gene coverage
# -----------------------------------------------------------------------------
cat("Step 3: Checking PAM50 gene coverage...\n")

available_genes <- pam50_genes_rat[pam50_genes_rat %in% rownames(pam50_data)]
missing_genes <- pam50_genes_human[!pam50_genes_rat %in% rownames(pam50_data)]

cat("  Present in data:", length(available_genes), "\n")
cat("  Missing from data:", length(missing_genes), "\n")

if (length(missing_genes) > 0) {
  cat("  Missing genes:", paste(head(missing_genes, 10), collapse = ", "))
  if (length(missing_genes) > 10) cat("...")
  cat("\n")
}

# Warn if too many genes missing
if (length(available_genes) < 40) {
  warning("Only ", length(available_genes), " of 50 PAM50 genes available. Results may be unreliable.")
}

# -----------------------------------------------------------------------------
# Step 4: Perform PAM50 classification
# -----------------------------------------------------------------------------
cat("Step 4: Running PAM50 classification...\n")

# Prepare expression matrix
pam50_tpm <- t(pam50_data[available_genes, ])
colnames(pam50_tpm) <- names(available_genes)

# FIX: Create minimal annotation data frame instead of using annot.nkis
# which is designed for human breast cancer samples
# The annotation only needs probe/gene identifiers for mapping
annot_df <- data.frame(
  probe = colnames(pam50_tpm),
  Gene.Symbol = colnames(pam50_tpm),
  EntrezGene.ID = NA,
  row.names = colnames(pam50_tpm)
)

PAM50_subtype <- molecular.subtyping(
  sbt.model = "pam50",
  data = pam50_tpm,
  annot = annot_df,
  do.mapping = FALSE  # Genes already mapped to human symbols
)

cat("  Subtype distribution:\n")
print(table(PAM50_subtype$subtype))

# -----------------------------------------------------------------------------
# Step 5: Save results
# -----------------------------------------------------------------------------
cat("Step 5: Saving results...\n")

# Save subtypes
subtype_df <- data.frame(
  Sample = names(PAM50_subtype$subtype),
  Subtype = PAM50_subtype$subtype,
  stringsAsFactors = FALSE
)
subtype_proba <- as.data.frame(PAM50_subtype$subtype.proba)
subtype_df <- cbind(subtype_df, subtype_proba)

write.csv(subtype_df, file.path(output_dir, "pam50_subtypes.csv"), row.names = FALSE)
cat("  Saved: pam50_subtypes.csv\n")

# Generate heatmap
centered_data <- t(scale(t(pam50_data[available_genes, ])))
annotation_col <- data.frame(
  Subtype = PAM50_subtype$subtype,
  row.names = colnames(centered_data)
)

pdf(file.path(output_dir, "pam50_heatmap.pdf"), width = 12, height = 15)
pheatmap(
  centered_data,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "PAM50 Gene Expression with Predicted Subtypes",
  fontsize_row = 10,
  fontsize_col = 10,
  cluster_cols = FALSE
)
dev.off()
cat("  Saved: pam50_heatmap.pdf\n")

# Save probability heatmap
pdf(file.path(output_dir, "pam50_probabilities_heatmap.pdf"), width = 12, height = 8)
pheatmap(
  t(PAM50_subtype$subtype.proba),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "PAM50 Subtype Probabilities",
  fontsize_row = 10,
  cluster_cols = FALSE,
  cluster_rows = FALSE
)
dev.off()
cat("  Saved: pam50_probabilities_heatmap.pdf\n")

cat("\n=== PAM50 complete ===\n")
