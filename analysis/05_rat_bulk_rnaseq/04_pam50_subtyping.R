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

# Set Arial as default font for all plots
library(showtext)
font_add("Arial", "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Arial.ttf")
showtext_auto()
theme_set(theme_bw(base_size = 14, base_family = "Arial"))

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
  "CDC20" = "Cdc20", "CDCA1" = "Nuf2", "CDC6" = "Cdc6",
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

# Load size-factor-normalized CPM data from DESeq2
cpm_file <- file.path(output_dir, "normalized_cpm.csv")
if (file.exists(cpm_file)) {
  pam50_data <- read.csv(cpm_file, row.names = 1)
} else {
  stop("Normalized CPM file not found. Run normalization first: ", cpm_file)
}

cat("  Expression matrix:", nrow(pam50_data), "genes x", ncol(pam50_data), "samples\n")

# -----------------------------------------------------------------------------
# Step 2b: Convert Ensembl IDs to gene symbols
# -----------------------------------------------------------------------------
cat("Step 2b: Converting Ensembl IDs to gene symbols...\n")

# Check if row names are Ensembl IDs
if (grepl("^ENSRNOG", rownames(pam50_data)[1])) {
  cat("  Detected Ensembl IDs, fetching gene symbols from biomaRt...\n")

  cache_file <- file.path(output_dir, "ensembl_to_symbol_cache.csv")
  if (file.exists(cache_file)) {
    cat("  Loading cached Ensembl→symbol mapping\n")
    gene_map <- read.csv(cache_file, stringsAsFactors = FALSE)
  } else {
    # Try main Ensembl, then mirrors
    gene_map <- NULL
    for (host in c("https://www.ensembl.org", "https://useast.ensembl.org", "https://asia.ensembl.org")) {
      gene_map <- tryCatch({
        mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = host)
        getBM(
          filters = "ensembl_gene_id",
          attributes = c("ensembl_gene_id", "external_gene_name"),
          values = rownames(pam50_data),
          mart = mart
        )
      }, error = function(e) {
        cat(sprintf("  BioMart mirror %s failed: %s\n", host, conditionMessage(e)))
        NULL
      })
      if (!is.null(gene_map) && nrow(gene_map) > 0) {
        cat(sprintf("  BioMart success via %s (%d mappings)\n", host, nrow(gene_map)))
        write.csv(gene_map, cache_file, row.names = FALSE)
        break
      }
    }
    if (is.null(gene_map) || nrow(gene_map) == 0) {
      stop("All BioMart mirrors failed. Retry later or provide a cached mapping at:\n  ", cache_file)
    }
  }

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
# Log2-transform CPM (genefu expects log2-scale data)
pam50_log2 <- log2(pam50_data[available_genes, ] + 1)

# Median-center per gene across samples (genefu PAM50 requirement:
# correlations to centroids assume centered data)
pam50_centered <- sweep(pam50_log2, 1, apply(pam50_log2, 1, median), "-")

# Transpose: genefu expects samples as rows, genes as columns
pam50_input <- t(pam50_centered)
colnames(pam50_input) <- names(available_genes)

# FIX: Create minimal annotation data frame instead of using annot.nkis
# which is designed for human breast cancer samples
# The annotation only needs probe/gene identifiers for mapping
annot_df <- data.frame(
  probe = colnames(pam50_input),
  Gene.Symbol = colnames(pam50_input),
  EntrezGene.ID = NA,
  row.names = colnames(pam50_input)
)

PAM50_subtype <- molecular.subtyping(
  sbt.model = "pam50",
  data = pam50_input,
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
# Fix R-mangled sample names (X102.FF.Tumor → 102-FF-Tumor)
colnames(centered_data) <- gsub("^X", "", gsub("\\.", "-", colnames(centered_data)))
annotation_col <- data.frame(
  Subtype = PAM50_subtype$subtype,
  row.names = colnames(centered_data)
)

annotation_colors <- list(
  Subtype = c(Basal = "#D55E00", Her2 = "#0072B2", LumA = "#009E73", LumB = "#E69F00", Normal = "#CC79A7")
)

pdf(file.path(output_dir, "pam50_heatmap.pdf"), width = 12, height = 15)
pheatmap(
  centered_data,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("#0072B2", "#F0E442", "#D55E00"))(101),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "PAM50 Gene Expression with Predicted Subtypes",
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 45,
  cluster_cols = FALSE,
  fontfamily = "Arial"
)
dev.off()
cat("  Saved: pam50_heatmap.pdf\n")

figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
png(file.path(figures_dir, "pam50_heatmap.png"), width = 12*300, height = 15*300, res = 300)
pheatmap(
  centered_data,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("#0072B2", "#F0E442", "#D55E00"))(101),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "PAM50 Gene Expression with Predicted Subtypes",
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 45,
  cluster_cols = FALSE,
  fontfamily = "Arial"
)
dev.off()
svg(file.path(figures_dir, "pam50_heatmap.svg"), width = 12, height = 15)
pheatmap(
  centered_data,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("#0072B2", "#F0E442", "#D55E00"))(101),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "PAM50 Gene Expression with Predicted Subtypes",
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 45,
  cluster_cols = FALSE,
  fontfamily = "Arial"
)
dev.off()  # SVG for vector assembly
cat("  Saved: pam50_heatmap.png\n")

# Save probability heatmap
pdf(file.path(output_dir, "pam50_probabilities_heatmap.pdf"), width = 12, height = 8)
pheatmap(
  t(PAM50_subtype$subtype.proba),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "PAM50 Subtype Probabilities",
  color = colorRampPalette(c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00"))(101),
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 45,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontfamily = "Arial"
)
dev.off()
cat("  Saved: pam50_probabilities_heatmap.pdf\n")

png(file.path(figures_dir, "pam50_probabilities_heatmap.png"), width = 12*300, height = 8*300, res = 300)
pheatmap(
  t(PAM50_subtype$subtype.proba),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "PAM50 Subtype Probabilities",
  color = colorRampPalette(c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00"))(101),
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 45,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontfamily = "Arial"
)
dev.off()
svg(file.path(figures_dir, "pam50_probabilities_heatmap.svg"), width = 12, height = 8)
pheatmap(
  t(PAM50_subtype$subtype.proba),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "PAM50 Subtype Probabilities",
  color = colorRampPalette(c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00"))(101),
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 45,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  fontfamily = "Arial"
)
dev.off()  # SVG for vector assembly
cat("  Saved: pam50_probabilities_heatmap.png\n")

cat("\n=== PAM50 complete ===\n")
