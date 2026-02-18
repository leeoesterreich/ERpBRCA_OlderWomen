# Code Integration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Integrate human scRNA-seq and rat bulk RNA-seq code into ERpBRCA_OlderWomen repo with biostatistical review.

**Architecture:** Sequential review-then-integrate approach. Review existing code for bugs and p-value issues, create missing GSVA/PROGENy/WCSEA code for scRNA-seq, convert all to numbered scripts matching existing repo pattern.

**Tech Stack:** R (Seurat, GSVA, progeny, indepthPathway, DESeq2, genefu), Bash (STAR, HTSeq), SLURM

---

## Phase 1: Rat Bulk RNA-seq Integration

### Task 1: Create directory structure and fix Alignment.sh

**Files:**
- Create: `analysis/05_rat_bulk_rnaseq/`
- Create: `analysis/05_rat_bulk_rnaseq/01_alignment.sh`

**Step 1: Create directory**

```bash
mkdir -p analysis/05_rat_bulk_rnaseq/outputs
mkdir -p analysis/05_rat_bulk_rnaseq/logs
```

**Step 2: Create fixed 01_alignment.sh**

Fix bugs:
- `input_dir1` → `input_dir` (line 49)
- Add `#SBATCH --mail-type=FAIL` and `#SBATCH --mail-user=`
- Parameterize hardcoded paths

```bash
#!/bin/bash
#SBATCH --job-name=rat_alignment
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Rat bulk RNA-seq alignment with STAR
# Inputs: Trimmed FASTQ files
# Outputs: Sorted BAM files

set -euo pipefail

module purge
module load gcc/8.2.0
module load star/2.7.9a

# Configurable paths
PROJECT_ROOT="${PROJECT_ROOT:-$(dirname $(dirname $(realpath $0)))}"
INPUT_DIR="${INPUT_DIR:-$PROJECT_ROOT/data/rat_bulk_rnaseq/trimmed}"
OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT_ROOT/analysis/05_rat_bulk_rnaseq/outputs/aligned}"
GENOME_DIR="${GENOME_DIR:-/path/to/Reference_genome/Rat}"
GENOME_FASTA="${GENOME_DIR}/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa"
GTF_FILE="${GENOME_DIR}/Rattus_norvegicus.mRatBN7.2.112.gtf"
THREADS=64

mkdir -p "$OUTPUT_DIR"

echo "=== STAR Alignment ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"

# Genome indexing (skip if index exists)
if [ ! -f "${GENOME_DIR}/SAindex" ]; then
    echo "Indexing reference genome..."
    STAR --runMode genomeGenerate \
         --genomeDir "$GENOME_DIR" \
         --genomeFastaFiles "$GENOME_FASTA" \
         --sjdbGTFfile "$GTF_FILE" \
         --runThreadN $THREADS
fi

# Align each sample
# FIX: Changed input_dir1 to INPUT_DIR
for read1 in "$INPUT_DIR"/*_R1_paired.fastq; do
    [ -f "$read1" ] || continue

    sample=$(basename "$read1" _R1_paired.fastq)
    read2="$INPUT_DIR/${sample}_R2_paired.fastq"
    output_prefix="$OUTPUT_DIR/${sample}_"

    if [ ! -f "$read2" ]; then
        echo "Warning: R2 not found for $sample, skipping"
        continue
    fi

    echo "Aligning: $sample"
    STAR --genomeDir "$GENOME_DIR" \
         --readFilesIn "$read1" "$read2" \
         --runThreadN $THREADS \
         --outFileNamePrefix "$output_prefix" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts

    echo "Completed: $sample"
done

echo "=== Alignment complete ==="
```

**Step 3: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/
git commit -m "feat(05_rat_bulk): add alignment script with bug fixes

Fixes:
- input_dir1 -> INPUT_DIR variable mismatch
- Added SBATCH mail directives
- Parameterized hardcoded paths"
```

---

### Task 2: Create fixed HTSeq script

**Files:**
- Create: `analysis/05_rat_bulk_rnaseq/02_htseq_count.sh`

**Step 1: Create 02_htseq_count.sh with strand flag**

```bash
#!/bin/bash
#SBATCH --job-name=rat_htseq
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# HTSeq read counting for rat bulk RNA-seq
# Inputs: Aligned BAM files
# Outputs: Gene count files

set -euo pipefail

module purge
module load htseq/0.13.5

PROJECT_ROOT="${PROJECT_ROOT:-$(dirname $(dirname $(realpath $0)))}"
INPUT_DIR="${INPUT_DIR:-$PROJECT_ROOT/analysis/05_rat_bulk_rnaseq/outputs/aligned}"
OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT_ROOT/analysis/05_rat_bulk_rnaseq/outputs/counts}"
GTF_FILE="${GTF_FILE:-/path/to/Rattus_norvegicus.mRatBN7.2.112.gtf}"

# FIX: Added strand specification - check library prep protocol
# Options: yes (forward), reverse, no (unstranded)
STRAND="${STRAND:-reverse}"

mkdir -p "$OUTPUT_DIR"

echo "=== HTSeq Count ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Strand: $STRAND"

for bam_file in "$INPUT_DIR"/*Aligned.sortedByCoord.out.bam; do
    [ -f "$bam_file" ] || continue

    sample=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)
    output_file="${OUTPUT_DIR}/${sample}.txt"

    echo "Counting: $sample"
    # FIX: Added -s flag for strand specification
    htseq-count -f bam -r pos -s "$STRAND" "$bam_file" "$GTF_FILE" > "$output_file"

    echo "Created: $output_file"
done

echo "=== HTSeq complete ==="
```

**Step 2: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/02_htseq_count.sh
git commit -m "feat(05_rat_bulk): add HTSeq counting with strand flag

FIX: Added -s flag for strand specification (default: reverse)"
```

---

### Task 3: Create fixed DESeq2 script

**Files:**
- Create: `analysis/05_rat_bulk_rnaseq/03_deseq2.R`

**Step 1: Create 03_deseq2.R with FDR filtering**

```r
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
```

**Step 2: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/03_deseq2.R
git commit -m "feat(05_rat_bulk): add DESeq2 script with FDR filtering

BIOSTATISTICAL FIX: Filter significant genes by padj (FDR) < 0.05
instead of raw p-value"
```

---

### Task 4: Create fixed PAM50 subtyping script

**Files:**
- Create: `analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R`

**Step 1: Create 04_pam50_subtyping.R with annotation fix**

```r
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
```

**Step 2: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R
git commit -m "feat(05_rat_bulk): add PAM50 subtyping with annotation fix

BIOSTATISTICAL FIX: Create minimal annotation instead of using
annot.nkis which is designed for human samples. Added gene
coverage check and warning when too many genes missing."
```

---

### Task 5: Create run scripts for rat bulk

**Files:**
- Create: `analysis/05_rat_bulk_rnaseq/run_analysis.sh`
- Create: `analysis/05_rat_bulk_rnaseq/run_analysis.sbatch`

**Step 1: Create run_analysis.sh**

```bash
#!/bin/bash
# Run rat bulk RNA-seq analysis pipeline
# Usage: bash run_analysis.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Rat Bulk RNA-seq Pipeline ==="
echo "Directory: $SCRIPT_DIR"

# Step 1: Alignment (submit as SLURM job)
echo "Step 1: Submitting alignment job..."
JOB1=$(sbatch --parsable 01_alignment.sh)
echo "  Job ID: $JOB1"

# Step 2: HTSeq counting (depends on alignment)
echo "Step 2: Submitting HTSeq job..."
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 02_htseq_count.sh)
echo "  Job ID: $JOB2"

# Step 3-4: R analysis (submit after counting)
echo "Step 3-4: Submitting R analysis job..."
sbatch --dependency=afterok:$JOB2 run_analysis.sbatch

echo "=== Jobs submitted ==="
```

**Step 2: Create run_analysis.sbatch**

```bash
#!/bin/bash
#SBATCH --job-name=rat_bulk_analysis
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/analysis_%j.out
#SBATCH --error=logs/analysis_%j.err

# Rat bulk RNA-seq R analysis
set -euo pipefail

module purge
module load r/4.2.0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p logs

echo "=== DESeq2 Analysis ==="
Rscript 03_deseq2.R

echo "=== PAM50 Subtyping ==="
Rscript 04_pam50_subtyping.R

echo "=== Analysis complete ==="
```

**Step 3: Commit**

```bash
chmod +x analysis/05_rat_bulk_rnaseq/run_analysis.sh
git add analysis/05_rat_bulk_rnaseq/run_analysis.sh analysis/05_rat_bulk_rnaseq/run_analysis.sbatch
git commit -m "feat(05_rat_bulk): add run scripts matching repo pattern"
```

---

## Phase 2: Human scRNA-seq Integration

### Task 6: Create directory and first scripts

**Files:**
- Create: `analysis/04_human_scrnaseq/`
- Create: `analysis/04_human_scrnaseq/01_load_subset_data.R`

**Step 1: Create directory**

```bash
mkdir -p analysis/04_human_scrnaseq/outputs
mkdir -p analysis/04_human_scrnaseq/logs
```

**Step 2: Create 01_load_subset_data.R**

```r
#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/01_load_subset_data.R
# Load and subset Wu et al. scRNA-seq data for age group analysis
#
# Inputs:
#   - data/human_scrnaseq/SeuratObj_GSE176078_ERpos_NewMeta_AfterQCSCT.rds
#   - data/human_scrnaseq/ClinicalData_Wu_scRNAseq_26p.txt
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_young_midage_elderly.rds

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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
data_dir <- file.path(project_root, "data/human_scrnaseq")

cat("=== Load and Subset scRNA-seq Data ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load clinical data
# -----------------------------------------------------------------------------
cat("Step 1: Loading clinical data...\n")

clinical_file <- file.path(data_dir, "ClinicalData_Wu_scRNAseq_26p.txt")
clinical_data <- fread(clinical_file, header = TRUE, stringsAsFactors = FALSE)
colnames(clinical_data) <- gsub(" ", "", colnames(clinical_data))

# Define age groups
clinical_data$AgeGroup <- ifelse(
  clinical_data$Age > 80, "Elderly",
  ifelse(clinical_data$Age <= 50, "Young", "MidAge")
)

cat("  Total patients:", nrow(clinical_data), "\n")
cat("  Age groups:\n")
print(table(clinical_data$AgeGroup))

# -----------------------------------------------------------------------------
# Step 2: Load Seurat object
# -----------------------------------------------------------------------------
cat("\nStep 2: Loading Seurat object...\n")

seurat_file <- file.path(data_dir, "SeuratObj_GSE176078_ERpos_NewMeta_AfterQCSCT.rds")
seurat_obj <- readRDS(seurat_file)

cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

# Run PCA and UMAP
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30, verbose = FALSE)

# -----------------------------------------------------------------------------
# Step 3: Subset to Young/MidAge/Elderly
# -----------------------------------------------------------------------------
cat("\nStep 3: Subsetting to age groups...\n")

# Sample IDs for each group (from Wu et al. clinical data)
young_ids <- c("CID3941", "CID4530N", "CID4535")
midage_ids <- c("CID4463", "CID4040", "CID4471", "CID4461")
elderly_ids <- c("CID3948", "CID4067", "CID4290A")
all_ids <- c(young_ids, midage_ids, elderly_ids)

seurat_subset <- subset(seurat_obj, subset = orig.ident %in% all_ids)

# Add age group labels
seurat_subset$AgeGroup <- case_when(
  seurat_subset$orig.ident %in% young_ids ~ "Young",
  seurat_subset$orig.ident %in% midage_ids ~ "MidAge",
  seurat_subset$orig.ident %in% elderly_ids ~ "Elderly",
  TRUE ~ NA_character_
)

cat("  Subset cells:", ncol(seurat_subset), "\n")
cat("  By age group:\n")
print(table(seurat_subset$AgeGroup))

# -----------------------------------------------------------------------------
# Step 4: Save output
# -----------------------------------------------------------------------------
cat("\nStep 4: Saving output...\n")

saveRDS(seurat_subset, file.path(output_dir, "seurat_young_midage_elderly.rds"))
saveRDS(clinical_data, file.path(output_dir, "clinical_data.rds"))

cat("\n=== Load complete ===\n")
```

**Step 3: Commit**

```bash
git add analysis/04_human_scrnaseq/
git commit -m "feat(04_human_scrnaseq): add data loading script"
```

---

### Task 7: Create Harmony integration script

**Files:**
- Create: `analysis/04_human_scrnaseq/02_harmony_integrate.R`

**Step 1: Create script**

```r
#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/02_harmony_integrate.R
# Harmony batch correction across patients
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_young_midage_elderly.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_harmony_integrated.rds
#   - analysis/04_human_scrnaseq/outputs/umap_by_patient.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

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

cat("=== Harmony Integration ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_young_midage_elderly.rds"))
cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Run Harmony
# -----------------------------------------------------------------------------
cat("\nStep 2: Running Harmony integration...\n")

seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "orig.ident",
  reduction = "pca",
  dims.use = 1:30,
  verbose = FALSE
)

# Run UMAP on Harmony reduction
seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.5)

cat("  Clusters found:", length(unique(seurat_obj$seurat_clusters)), "\n")

# -----------------------------------------------------------------------------
# Step 3: Visualize
# -----------------------------------------------------------------------------
cat("\nStep 3: Generating visualizations...\n")

p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
  ggtitle("By Patient")

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "AgeGroup") +
  ggtitle("By Age Group")

p3 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
  ggtitle("Clusters")

pdf(file.path(output_dir, "umap_harmony.pdf"), width = 15, height = 5)
print(p1 + p2 + p3)
dev.off()

# -----------------------------------------------------------------------------
# Step 4: Save
# -----------------------------------------------------------------------------
cat("\nStep 4: Saving...\n")

saveRDS(seurat_obj, file.path(output_dir, "seurat_harmony_integrated.rds"))

cat("\n=== Harmony complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/04_human_scrnaseq/02_harmony_integrate.R
git commit -m "feat(04_human_scrnaseq): add Harmony integration"
```

---

### Task 8: Create cell type annotation script

**Files:**
- Create: `analysis/04_human_scrnaseq/03_cell_type_annotation.R`

**Step 1: Create script (adapting from original with bug fixes)**

```r
#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/03_cell_type_annotation.R
# Cell type annotation using marker genes
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_harmony_integrated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#   - analysis/04_human_scrnaseq/outputs/umap_celltypes.pdf
#   - analysis/04_human_scrnaseq/outputs/metadata_annotated.txt

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(patchwork)
})

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

cat("=== Cell Type Annotation ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_harmony_integrated.rds"))
cat("  Cells:", ncol(seurat_obj), "\n")

# Check if cell type annotations exist from original data
if ("CellTypeMajor" %in% colnames(seurat_obj@meta.data)) {
  cat("  Using existing cell type annotations\n")

  # Clean up cell type names (remove special characters)
  seurat_obj$CellTypeMajor <- gsub("-", "", seurat_obj$CellTypeMajor)

  # Create combined annotation
  # FIX: Ensure CellTypeMinor exists before using
  if ("CellTypeMinor" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$CellTypeAnnot <- ifelse(
      seurat_obj$CellTypeMajor %in% c("Myeloid", "Tcells"),
      seurat_obj$CellTypeMinor,
      seurat_obj$CellTypeMajor
    )
  } else {
    seurat_obj$CellTypeAnnot <- seurat_obj$CellTypeMajor
  }

  # Clean annotation names
  seurat_obj$CellTypeAnnot <- gsub("-|_|\\+", "", seurat_obj$CellTypeAnnot)

} else {
  cat("  No existing annotations, using cluster-based annotation\n")
  # Placeholder - would use marker-based annotation here
  seurat_obj$CellTypeAnnot <- paste0("Cluster_", seurat_obj$seurat_clusters)
}

cat("  Cell types:\n")
print(table(seurat_obj$CellTypeAnnot))

# -----------------------------------------------------------------------------
# Step 2: Set identity and visualize
# -----------------------------------------------------------------------------
cat("\nStep 2: Visualizing...\n")

Idents(seurat_obj) <- seurat_obj$CellTypeAnnot

p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4) +
  ggtitle("Cell Types") +
  theme(legend.position = "right")

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "AgeGroup") +
  ggtitle("Age Groups")

pdf(file.path(output_dir, "umap_celltypes.pdf"), width = 14, height = 6)
print(p1 + p2)
dev.off()

# Feature plots for key markers
markers <- c("EPCAM", "KRT19", "CD68", "CD3D", "MS4A1", "PECAM1")
markers_present <- markers[markers %in% rownames(seurat_obj)]

if (length(markers_present) > 0) {
  pdf(file.path(output_dir, "feature_markers.pdf"), width = 12, height = 8)
  print(FeaturePlot(seurat_obj, features = markers_present, ncol = 3))
  dev.off()
}

# -----------------------------------------------------------------------------
# Step 3: Save metadata
# -----------------------------------------------------------------------------
cat("\nStep 3: Saving...\n")

# Save annotated object
saveRDS(seurat_obj, file.path(output_dir, "seurat_annotated.rds"))

# Save metadata
metadata <- seurat_obj@meta.data %>%
  tibble::rownames_to_column("CellID") %>%
  select(CellID, orig.ident, AgeGroup, CellTypeAnnot, seurat_clusters)

fwrite(metadata, file.path(output_dir, "metadata_annotated.txt"),
       sep = "\t", quote = FALSE)

cat("\n=== Annotation complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/04_human_scrnaseq/03_cell_type_annotation.R
git commit -m "feat(04_human_scrnaseq): add cell type annotation

BUG FIX: Check CellTypeMinor exists before using"
```

---

### Task 9: Create cell fractions script with FDR correction

**Files:**
- Create: `analysis/04_human_scrnaseq/04_cell_fractions.R`

**Step 1: Create script with p-value correction**

```r
#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/04_cell_fractions.R
# Cell type fraction analysis across age groups
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#   - analysis/04_human_scrnaseq/outputs/clinical_data.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/cell_fractions.csv
#   - analysis/04_human_scrnaseq/outputs/cell_fraction_stats.csv
#   - analysis/04_human_scrnaseq/outputs/fraction_boxplot.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
})

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

cat("=== Cell Fraction Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))
clinical_data <- readRDS(file.path(output_dir, "clinical_data.rds"))

cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Calculate cell type fractions per patient
# -----------------------------------------------------------------------------
cat("\nStep 2: Calculating cell fractions...\n")

# Count cells per patient per cell type
cell_counts <- table(seurat_obj$CellTypeAnnot, seurat_obj$orig.ident) %>%
  as.data.frame.matrix()

# Calculate proportions
cell_fractions <- sweep(cell_counts, 2, colSums(cell_counts), "/") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("CaseID")

# Merge with clinical data
fraction_data <- inner_join(
  clinical_data[, c("CaseID", "Age", "AgeGroup")],
  cell_fractions,
  by = "CaseID"
)

write.csv(fraction_data, file.path(output_dir, "cell_fractions.csv"), row.names = FALSE)
cat("  Saved: cell_fractions.csv\n")

# -----------------------------------------------------------------------------
# Step 3: Statistical tests with FDR correction
# -----------------------------------------------------------------------------
cat("\nStep 3: Running statistical tests...\n")

cell_types <- colnames(cell_fractions)[-1]  # Exclude CaseID

# Compare Young vs Elderly for each cell type
# FIX: Apply FDR correction across all comparisons
stat_results <- lapply(cell_types, function(ct) {
  young_vals <- fraction_data %>% filter(AgeGroup == "Young") %>% pull(!!sym(ct))
  elderly_vals <- fraction_data %>% filter(AgeGroup == "Elderly") %>% pull(!!sym(ct))

  # Wilcoxon test (non-parametric, appropriate for small n)
  if (length(young_vals) >= 2 && length(elderly_vals) >= 2) {
    test_result <- wilcox.test(young_vals, elderly_vals, exact = FALSE)
    data.frame(
      CellType = ct,
      Young_mean = mean(young_vals),
      Elderly_mean = mean(elderly_vals),
      pvalue = test_result$p.value
    )
  } else {
    data.frame(
      CellType = ct,
      Young_mean = mean(young_vals),
      Elderly_mean = mean(elderly_vals),
      pvalue = NA
    )
  }
}) %>% bind_rows()

# BIOSTATISTICAL FIX: Apply BH-FDR correction
stat_results$padj <- p.adjust(stat_results$pvalue, method = "BH")
stat_results$significant <- stat_results$padj < 0.05

write.csv(stat_results, file.path(output_dir, "cell_fraction_stats.csv"), row.names = FALSE)

cat("  Cell types tested:", nrow(stat_results), "\n")
cat("  Significant (FDR < 0.05):", sum(stat_results$significant, na.rm = TRUE), "\n")

# -----------------------------------------------------------------------------
# Step 4: Visualization
# -----------------------------------------------------------------------------
cat("\nStep 4: Generating plots...\n")

# Reshape for plotting
plot_data <- fraction_data %>%
  pivot_longer(
    cols = all_of(cell_types),
    names_to = "CellType",
    values_to = "Fraction"
  ) %>%
  filter(AgeGroup != "MidAge")  # Compare Young vs Elderly

p <- ggplot(plot_data, aes(x = CellType, y = Fraction, fill = AgeGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = AgeGroup), position = position_jitterdodge(jitter.width = 0.1)) +
  scale_fill_manual(values = c("Young" = "#4DAF4A", "Elderly" = "#E41A1C")) +
  scale_color_manual(values = c("Young" = "#4DAF4A", "Elderly" = "#E41A1C")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(x = "", y = "Cell Type Fraction", title = "Cell Type Proportions by Age Group")

pdf(file.path(output_dir, "fraction_boxplot.pdf"), width = 12, height = 6)
print(p)
dev.off()

cat("\n=== Cell fraction analysis complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/04_human_scrnaseq/04_cell_fractions.R
git commit -m "feat(04_human_scrnaseq): add cell fraction analysis

BIOSTATISTICAL FIX: Added BH-FDR correction for multiple cell type
comparisons (Young vs Elderly)"
```

---

### Task 10: Create gene expression violin script

**Files:**
- Create: `analysis/04_human_scrnaseq/05_gene_expression_violin.R`

**Step 1: Create script**

```r
#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/05_gene_expression_violin.R
# Gene expression violin plots by age group and cell type
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/violin_*.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

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

cat("=== Gene Expression Violin Plots ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load and prepare data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))

# Normalize if not already done
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

# Order age groups
seurat_obj$AgeGroup <- factor(seurat_obj$AgeGroup, levels = c("Young", "MidAge", "Elderly"))

cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Define gene sets
# -----------------------------------------------------------------------------
cat("\nStep 2: Defining gene sets...\n")

immune_checkpoint <- c("CTLA4", "PDCD1", "CD274", "LAG3", "HAVCR2", "PDCD1LG2")
macrophage_markers <- c("CCL2", "CCL3", "CCL4", "TNF", "TGFB1", "CD163", "MRC1")
cytokines <- c("IL1A", "IL1B", "IL6", "CXCL8", "CXCL10", "CCL5")

# Filter to genes present in data
filter_genes <- function(genes) {
  present <- genes[genes %in% rownames(seurat_obj)]
  if (length(present) < length(genes)) {
    missing <- setdiff(genes, present)
    cat("  Missing:", paste(missing, collapse = ", "), "\n")
  }
  present
}

immune_checkpoint <- filter_genes(immune_checkpoint)
macrophage_markers <- filter_genes(macrophage_markers)
cytokines <- filter_genes(cytokines)

# -----------------------------------------------------------------------------
# Step 3: Generate violin plots by cell type
# -----------------------------------------------------------------------------
cat("\nStep 3: Generating plots...\n")

cell_types <- unique(seurat_obj$CellTypeAnnot)

for (ct in cell_types) {
  cat("  Processing:", ct, "\n")

  # Subset to cell type
  seurat_subset <- subset(seurat_obj, CellTypeAnnot == ct)

  if (ncol(seurat_subset) < 10) {
    cat("    Skipping (too few cells)\n")
    next
  }

  # Immune checkpoint genes
  if (length(immune_checkpoint) > 0) {
    p <- VlnPlot(
      seurat_subset,
      features = immune_checkpoint,
      group.by = "AgeGroup",
      pt.size = 0
    ) + plot_annotation(title = paste(ct, "- Immune Checkpoint"))

    pdf(file.path(output_dir, paste0("violin_checkpoint_", ct, ".pdf")), width = 12, height = 8)
    print(p)
    dev.off()
  }

  # Macrophage markers (only for macrophage/monocyte)
  if (grepl("Macro|Mono", ct, ignore.case = TRUE) && length(macrophage_markers) > 0) {
    p <- VlnPlot(
      seurat_subset,
      features = macrophage_markers,
      group.by = "AgeGroup",
      pt.size = 0
    ) + plot_annotation(title = paste(ct, "- Macrophage Markers"))

    pdf(file.path(output_dir, paste0("violin_macrophage_", ct, ".pdf")), width = 12, height = 8)
    print(p)
    dev.off()
  }
}

cat("\n=== Violin plots complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/04_human_scrnaseq/05_gene_expression_violin.R
git commit -m "feat(04_human_scrnaseq): add gene expression violin plots"
```

---

### Task 11: Create scRNA GSVA script (NEW - adapted from bulk)

**Files:**
- Create: `analysis/04_human_scrnaseq/06_run_gsva.R`

**Step 1: Create script adapted for single-cell (pseudo-bulk approach)**

```r
#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/06_run_gsva.R
# GSVA on single-cell data using pseudo-bulk approach
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/gsva_pseudobulk.rds
#   - analysis/04_human_scrnaseq/outputs/gsva_heatmap.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(GSVA)
  library(msigdbr)
  library(dplyr)
  library(data.table)
  library(pheatmap)
  library(tibble)
})

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

cat("=== scRNA-seq GSVA (Pseudo-bulk) ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))
cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Create pseudo-bulk per patient
# -----------------------------------------------------------------------------
cat("\nStep 2: Creating pseudo-bulk profiles...\n")

# Aggregate counts per patient
DefaultAssay(seurat_obj) <- "RNA"

# Get raw counts
counts <- GetAssayData(seurat_obj, slot = "counts")

# Aggregate by patient (orig.ident)
patients <- unique(seurat_obj$orig.ident)
pseudobulk <- sapply(patients, function(pt) {
  cells <- colnames(seurat_obj)[seurat_obj$orig.ident == pt]
  rowSums(counts[, cells, drop = FALSE])
})

# Normalize (CPM + log)
pseudobulk_cpm <- sweep(pseudobulk, 2, colSums(pseudobulk), "/") * 1e6
pseudobulk_log <- log2(pseudobulk_cpm + 1)

cat("  Pseudo-bulk matrix:", nrow(pseudobulk_log), "genes x", ncol(pseudobulk_log), "samples\n")

# -----------------------------------------------------------------------------
# Step 3: Load gene sets (same as bulk analysis)
# -----------------------------------------------------------------------------
cat("\nStep 3: Loading gene sets...\n")

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

# GO BP estrogen pathways
gobp_sets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gobp_list <- split(gobp_sets$gene_symbol, gobp_sets$gs_name)
gobp_estrogen <- gobp_list[c(
  "GOBP_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_CELLULAR_RESPONSE_TO_ESTROGEN_STIMULUS"
)]

# Combine
estrogen_pathways <- c(hallmark_estrogen, reactome_estrogen, gobp_estrogen)
estrogen_pathways <- estrogen_pathways[!sapply(estrogen_pathways, is.null)]
cat("  Gene sets:", length(estrogen_pathways), "\n")

# -----------------------------------------------------------------------------
# Step 4: Run GSVA
# -----------------------------------------------------------------------------
cat("\nStep 4: Running GSVA...\n")

gsva_result <- gsva(
  gsvaParam(
    as.matrix(pseudobulk_log),
    estrogen_pathways,
    kcdf = "Gaussian",
    maxDiff = TRUE
  )
)

cat("  Result:", nrow(gsva_result), "pathways x", ncol(gsva_result), "samples\n")

# -----------------------------------------------------------------------------
# Step 5: Add age group annotation and visualize
# -----------------------------------------------------------------------------
cat("\nStep 5: Visualizing...\n")

# Get age group per patient
age_groups <- seurat_obj@meta.data %>%
  select(orig.ident, AgeGroup) %>%
  distinct()
rownames(age_groups) <- age_groups$orig.ident

# Annotation for heatmap
ann_col <- data.frame(
  AgeGroup = age_groups[colnames(gsva_result), "AgeGroup"],
  row.names = colnames(gsva_result)
)

ann_colors <- list(
  AgeGroup = c(Young = "#4DAF4A", MidAge = "#377EB8", Elderly = "#E41A1C")
)

pdf(file.path(output_dir, "gsva_heatmap.pdf"), width = 10, height = 6)
pheatmap(
  gsva_result,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  main = "GSVA Estrogen Pathway Scores (Pseudo-bulk)",
  scale = "row",
  show_colnames = TRUE
)
dev.off()

# -----------------------------------------------------------------------------
# Step 6: Save
# -----------------------------------------------------------------------------
cat("\nStep 6: Saving...\n")

saveRDS(gsva_result, file.path(output_dir, "gsva_pseudobulk.rds"))
saveRDS(pseudobulk_log, file.path(output_dir, "pseudobulk_log2cpm.rds"))

cat("\n=== scRNA GSVA complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/04_human_scrnaseq/06_run_gsva.R
git commit -m "feat(04_human_scrnaseq): add GSVA for scRNA-seq

NEW: Creates pseudo-bulk profiles per patient then runs GSVA.
Adapted from bulk RNA-seq GSVA script (01_human_bulk_rnaseq/02_run_gsva.R)"
```

---

### Task 12: Create scRNA PROGENy script (NEW)

**Files:**
- Create: `analysis/04_human_scrnaseq/07_run_progeny.R`

**Step 1: Create script**

```r
#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/07_run_progeny.R
# PROGENy pathway activity on single-cell data
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/progeny_activity.rds
#   - analysis/04_human_scrnaseq/outputs/progeny_heatmap.pdf

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(progeny)
  library(dplyr)
  library(pheatmap)
  library(tibble)
})

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

cat("=== scRNA-seq PROGENy ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))

# Use normalized data
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

cat("  Cells:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Create pseudo-bulk per patient
# -----------------------------------------------------------------------------
cat("\nStep 2: Creating pseudo-bulk profiles...\n")

# Use the log-normalized data
norm_data <- GetAssayData(seurat_obj, slot = "data")

# Aggregate by patient (mean expression)
patients <- unique(seurat_obj$orig.ident)
pseudobulk <- sapply(patients, function(pt) {
  cells <- colnames(seurat_obj)[seurat_obj$orig.ident == pt]
  rowMeans(norm_data[, cells, drop = FALSE])
})

cat("  Pseudo-bulk matrix:", nrow(pseudobulk), "genes x", ncol(pseudobulk), "samples\n")

# -----------------------------------------------------------------------------
# Step 3: Run PROGENy
# -----------------------------------------------------------------------------
cat("\nStep 3: Running PROGENy...\n")

pathway_activity <- progeny(
  as.matrix(pseudobulk),
  scale = FALSE,
  organism = "Human",
  top = 100,
  perm = 1000
)

cat("  Result:", nrow(pathway_activity), "samples x", ncol(pathway_activity), "pathways\n")

# -----------------------------------------------------------------------------
# Step 4: Visualize
# -----------------------------------------------------------------------------
cat("\nStep 4: Visualizing...\n")

# Get age group per patient
age_groups <- seurat_obj@meta.data %>%
  select(orig.ident, AgeGroup) %>%
  distinct()
rownames(age_groups) <- age_groups$orig.ident

# Annotation
ann_row <- data.frame(
  AgeGroup = age_groups[rownames(pathway_activity), "AgeGroup"],
  row.names = rownames(pathway_activity)
)

ann_colors <- list(
  AgeGroup = c(Young = "#4DAF4A", MidAge = "#377EB8", Elderly = "#E41A1C")
)

pdf(file.path(output_dir, "progeny_heatmap.pdf"), width = 10, height = 8)
pheatmap(
  pathway_activity,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  main = "PROGENy Pathway Activity (Pseudo-bulk)",
  scale = "column",
  cluster_cols = FALSE
)
dev.off()

# Focus on Estrogen pathway
if ("Estrogen" %in% colnames(pathway_activity)) {
  estrogen_activity <- data.frame(
    Patient = rownames(pathway_activity),
    Estrogen = pathway_activity[, "Estrogen"],
    AgeGroup = ann_row$AgeGroup
  )
  write.csv(estrogen_activity, file.path(output_dir, "progeny_estrogen_activity.csv"),
            row.names = FALSE)
}

# -----------------------------------------------------------------------------
# Step 5: Save
# -----------------------------------------------------------------------------
cat("\nStep 5: Saving...\n")

saveRDS(pathway_activity, file.path(output_dir, "progeny_activity.rds"))

cat("\n=== scRNA PROGENy complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/04_human_scrnaseq/07_run_progeny.R
git commit -m "feat(04_human_scrnaseq): add PROGENy for scRNA-seq

NEW: Creates pseudo-bulk profiles per patient then runs PROGENy.
Adapted from bulk RNA-seq PROGENy script."
```

---

### Task 13: Create WCSEA script (NEW - from indepthPathway)

**Files:**
- Create: `analysis/04_human_scrnaseq/08_run_wcsea.R`

**Step 1: Create script based on indepthPathway documentation**

```r
#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/08_run_wcsea.R
# WCSEA pathway analysis for single-cell RNA-seq
# Using indepthPathway package (https://github.com/wangxlab/indepthPathway)
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/wcsea_results.rds
#   - analysis/04_human_scrnaseq/outputs/wcsea_enrichment.csv

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

# Check if indepthPathway is installed
if (!requireNamespace("indepthPathway", quietly = TRUE)) {
  cat("Installing indepthPathway from GitHub...\n")
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("wangxlab/indepthPathway")
}
library(indepthPathway)

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

cat("=== WCSEA Pathway Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_obj <- readRDS(file.path(output_dir, "seurat_annotated.rds"))

# Get normalized expression matrix
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
expr_matrix <- GetAssayData(seurat_obj, slot = "data")

cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# Step 2: Identify marker genes for comparison
# -----------------------------------------------------------------------------
cat("\nStep 2: Finding markers (Young vs Elderly)...\n")

# Subset to Young and Elderly only
seurat_subset <- subset(seurat_obj, AgeGroup %in% c("Young", "Elderly"))
Idents(seurat_subset) <- seurat_subset$AgeGroup

# Find differentially expressed genes
markers <- FindMarkers(
  seurat_subset,
  ident.1 = "Elderly",
  ident.2 = "Young",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# Apply FDR correction (already done by FindMarkers via p_val_adj)
markers_sig <- markers %>%
  filter(p_val_adj < 0.05) %>%
  tibble::rownames_to_column("gene")

cat("  Significant markers:", nrow(markers_sig), "\n")
cat("    Upregulated in Elderly:", sum(markers_sig$avg_log2FC > 0), "\n")
cat("    Downregulated in Elderly:", sum(markers_sig$avg_log2FC < 0), "\n")

# -----------------------------------------------------------------------------
# Step 3: Run WCSEA
# -----------------------------------------------------------------------------
cat("\nStep 3: Running WCSEA...\n")

# Prepare gene list with fold changes
gene_fc <- setNames(markers_sig$avg_log2FC, markers_sig$gene)

# Run WCSEA (Weighted Concept Signature Enrichment Analysis)
# This is designed to handle noisy single-cell data
tryCatch({
  wcsea_result <- WCSEA(
    gene_list = gene_fc,
    species = "human",
    pathway_db = "GO_BP",  # Gene Ontology Biological Process
    min_genes = 10,
    max_genes = 500
  )

  # Extract results
  wcsea_df <- wcsea_result$enrichment_result

  # FDR correction
  wcsea_df$padj <- p.adjust(wcsea_df$pvalue, method = "BH")
  wcsea_sig <- wcsea_df %>%
    filter(padj < 0.05) %>%
    arrange(padj)

  cat("  Enriched pathways (FDR < 0.05):", nrow(wcsea_sig), "\n")

  # Save results
  saveRDS(wcsea_result, file.path(output_dir, "wcsea_results.rds"))
  write.csv(wcsea_sig, file.path(output_dir, "wcsea_enrichment.csv"), row.names = FALSE)

}, error = function(e) {
  cat("  Error running WCSEA:", conditionMessage(e), "\n")
  cat("  Saving marker genes for manual WCSEA analysis\n")

  # Save markers for manual analysis
  write.csv(markers_sig, file.path(output_dir, "wcsea_input_markers.csv"), row.names = FALSE)
})

# -----------------------------------------------------------------------------
# Step 4: Also check estrogen-related pathways specifically
# -----------------------------------------------------------------------------
cat("\nStep 4: Checking estrogen pathway genes...\n")

estrogen_genes <- c("ESR1", "ESR2", "PGR", "GREB1", "TFF1", "AREG", "XBP1")
estrogen_in_markers <- markers_sig %>%
  filter(gene %in% estrogen_genes)

if (nrow(estrogen_in_markers) > 0) {
  cat("  Estrogen-related genes in markers:\n")
  print(estrogen_in_markers[, c("gene", "avg_log2FC", "p_val_adj")])
} else {
  cat("  No canonical estrogen genes in significant markers\n")
}

cat("\n=== WCSEA complete ===\n")
```

**Step 2: Commit**

```bash
git add analysis/04_human_scrnaseq/08_run_wcsea.R
git commit -m "feat(04_human_scrnaseq): add WCSEA pathway analysis

NEW: Implements WCSEA from indepthPathway package for single-cell
pathway enrichment. Uses differential expression between age groups
as input. FDR correction applied to enrichment results."
```

---

### Task 14: Create CellPhoneDB prep script (with bug fixes)

**Files:**
- Create: `analysis/04_human_scrnaseq/09_cellphonedb_prep.R`

**Step 1: Create script with bug fixes from original**

```r
#!/usr/bin/env Rscript
# analysis/04_human_scrnaseq/09_cellphonedb_prep.R
# Prepare input files for CellPhoneDB analysis
#
# Inputs:
#   - analysis/04_human_scrnaseq/outputs/seurat_annotated.rds
#
# Outputs:
#   - analysis/04_human_scrnaseq/outputs/cellphonedb/counts_*.txt
#   - analysis/04_human_scrnaseq/outputs/cellphonedb/metadata_*.txt

set.seed(12345)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "outputs/cellphonedb")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== CellPhoneDB Input Preparation ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load data
# -----------------------------------------------------------------------------
cat("Step 1: Loading data...\n")

seurat_file <- file.path(dirname(output_dir), "seurat_annotated.rds")
seurat_obj <- readRDS(seurat_file)

cat("  Cells:", ncol(seurat_obj), "\n")

# FIX: Define these subsets properly instead of using undefined variables
seurat_young <- subset(seurat_obj, AgeGroup == "Young")
seurat_elderly <- subset(seurat_obj, AgeGroup == "Elderly")
seurat_whole <- seurat_obj

# -----------------------------------------------------------------------------
# Step 2: Function to prepare CellPhoneDB inputs
# -----------------------------------------------------------------------------
prepare_cellphonedb <- function(seurat_subset, name) {
  cat(sprintf("\nProcessing %s (%d cells)...\n", name, ncol(seurat_subset)))

  # Get count data
  counts <- GetAssayData(seurat_subset, slot = "counts")

  # Filter low-count genes (>20 total counts across cells)
  counts_filtered <- counts[rowSums(counts) > 20, ]

  # Filter low-count cells (>1000 total counts)
  counts_filtered <- counts_filtered[, colSums(counts_filtered) > 1000]

  # FIX: Replace hyphens in cell IDs (CellPhoneDB requirement)
  colnames(counts_filtered) <- gsub("-", "_", colnames(counts_filtered))

  # Save counts
  count_file <- file.path(output_dir, sprintf("counts_%s.txt", name))
  fwrite(
    as.data.frame(counts_filtered) %>% tibble::rownames_to_column("Gene"),
    count_file,
    sep = "\t",
    quote = FALSE
  )
  cat(sprintf("  Saved: %s (%d genes x %d cells)\n",
              basename(count_file), nrow(counts_filtered), ncol(counts_filtered)))

  # Prepare metadata
  cells_keep <- colnames(counts_filtered)
  # FIX: Match cell names after hyphen replacement
  original_cells <- gsub("_", "-", cells_keep)

  metadata <- data.frame(
    Cell = cells_keep,
    cell_type = gsub("-|_|\\+", "", seurat_subset$CellTypeAnnot[match(original_cells, colnames(seurat_subset))])
  )

  # FIX: Verify cell IDs match between counts and metadata
  if (!all(metadata$Cell == colnames(counts_filtered))) {
    warning("Cell ID mismatch between counts and metadata!")
  }

  # Save metadata
  meta_file <- file.path(output_dir, sprintf("metadata_%s.txt", name))
  fwrite(metadata, meta_file, sep = "\t", quote = FALSE)
  cat(sprintf("  Saved: %s\n", basename(meta_file)))

  return(list(
    n_genes = nrow(counts_filtered),
    n_cells = ncol(counts_filtered)
  ))
}

# -----------------------------------------------------------------------------
# Step 3: Prepare files for each group
# -----------------------------------------------------------------------------
cat("\nStep 3: Preparing CellPhoneDB inputs...\n")

stats_young <- prepare_cellphonedb(seurat_young, "Young")
stats_elderly <- prepare_cellphonedb(seurat_elderly, "Elderly")
stats_whole <- prepare_cellphonedb(seurat_whole, "Whole")

# -----------------------------------------------------------------------------
# Step 4: Summary
# -----------------------------------------------------------------------------
cat("\n=== Summary ===\n")
cat(sprintf("Young: %d genes x %d cells\n", stats_young$n_genes, stats_young$n_cells))
cat(sprintf("Elderly: %d genes x %d cells\n", stats_elderly$n_genes, stats_elderly$n_cells))
cat(sprintf("Whole: %d genes x %d cells\n", stats_whole$n_genes, stats_whole$n_cells))

cat("\n=== CellPhoneDB prep complete ===\n")
cat("Run CellPhoneDB with:\n")
cat("  cellphonedb method statistical_analysis metadata_<group>.txt counts_<group>.txt\n")
```

**Step 2: Commit**

```bash
git add analysis/04_human_scrnaseq/09_cellphonedb_prep.R
git commit -m "feat(04_human_scrnaseq): add CellPhoneDB prep

BUG FIXES:
- Define seurat_young/seurat_elderly properly (were undefined)
- Fix typo: CountData_YoungElderlyroWhole -> counts_filtered
- Verify cell ID consistency between counts and metadata"
```

---

### Task 15: Create run scripts for human scRNA

**Files:**
- Create: `analysis/04_human_scrnaseq/run_analysis.sh`
- Create: `analysis/04_human_scrnaseq/run_analysis.sbatch`

**Step 1: Create run_analysis.sh**

```bash
#!/bin/bash
# Run human scRNA-seq analysis pipeline
# Usage: bash run_analysis.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Human scRNA-seq Pipeline ==="
echo "Directory: $SCRIPT_DIR"

# Submit SLURM job
sbatch run_analysis.sbatch

echo "=== Job submitted ==="
```

**Step 2: Create run_analysis.sbatch**

```bash
#!/bin/bash
#SBATCH --job-name=human_scrnaseq
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/analysis_%j.out
#SBATCH --error=logs/analysis_%j.err

# Human scRNA-seq analysis pipeline
set -euo pipefail

module purge
module load r/4.2.0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p logs outputs

echo "=== Step 1: Load and subset data ==="
Rscript 01_load_subset_data.R

echo "=== Step 2: Harmony integration ==="
Rscript 02_harmony_integrate.R

echo "=== Step 3: Cell type annotation ==="
Rscript 03_cell_type_annotation.R

echo "=== Step 4: Cell fractions ==="
Rscript 04_cell_fractions.R

echo "=== Step 5: Gene expression ==="
Rscript 05_gene_expression_violin.R

echo "=== Step 6: GSVA ==="
Rscript 06_run_gsva.R

echo "=== Step 7: PROGENy ==="
Rscript 07_run_progeny.R

echo "=== Step 8: WCSEA ==="
Rscript 08_run_wcsea.R

echo "=== Step 9: CellPhoneDB prep ==="
Rscript 09_cellphonedb_prep.R

echo "=== Analysis complete ==="
```

**Step 3: Commit**

```bash
chmod +x analysis/04_human_scrnaseq/run_analysis.sh
git add analysis/04_human_scrnaseq/run_analysis.sh analysis/04_human_scrnaseq/run_analysis.sbatch
git commit -m "feat(04_human_scrnaseq): add run scripts matching repo pattern"
```

---

## Phase 3: Documentation

### Task 16: Update README and environment

**Files:**
- Modify: `README.md`
- Modify: `environment.yml`

**Step 1: Update README.md**

Add to repository structure section:

```markdown
├── analysis/
│   ├── 01_human_bulk_rnaseq/    # Human ER+ breast cancer RNA-seq
│   ├── 02_rat_snrnaseq/         # Rat mammary tumor snRNA-seq
│   ├── 03_comparison/           # Cross-species comparison
│   ├── 03_rat_wes/              # Rat whole exome sequencing
│   ├── 04_human_scrnaseq/       # Human ER+ scRNA-seq (Wu et al.)
│   └── 05_rat_bulk_rnaseq/      # Rat mammary tumor bulk RNA-seq
```

Add to Data Availability table:

```markdown
| Human scRNA-seq | [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) | Wu et al. 2021 |
```

Add to Methods Summary:

```markdown
### Human scRNA-seq (Wu et al. 2021)
- Seurat v5 + Harmony batch correction
- Pseudo-bulk GSVA and PROGENy analysis
- WCSEA pathway enrichment (indepthPathway)
- CellPhoneDB interaction analysis

### Rat Bulk RNA-seq
- STAR alignment to mRatBN7.2
- DESeq2 differential expression with BH-FDR correction
- PAM50 molecular subtyping (human orthologs)
```

**Step 2: Update environment.yml**

Add new dependencies:

```yaml
dependencies:
  # ... existing ...
  - r-genefu
  - bioconda::bioconductor-progeny
  - bioconda::bioconductor-gsva
  # indepthPathway installed from GitHub in script
```

**Step 3: Commit**

```bash
git add README.md environment.yml
git commit -m "docs: update README and environment for new analyses

Added:
- 04_human_scrnaseq (Wu et al. scRNA-seq)
- 05_rat_bulk_rnaseq (Rat tumor bulk RNA-seq)
- genefu dependency for PAM50 subtyping"
```

---

### Task 17: Create issues log

**Files:**
- Create: `docs/plans/2026-02-17-integration-issues-log.md`

**Step 1: Create issues log documenting all fixes**

```markdown
# Code Integration Issues Log

## Summary

| Category | Issues Found | Fixed |
|----------|-------------|-------|
| P-value correction | 3 | 3 |
| Code bugs | 6 | 6 |
| Missing code | 3 | 3 |
| HPC policy | 2 | 2 |

## Issues Fixed

### Human scRNA-seq (GitHub)

| Issue | File | Fix |
|-------|------|-----|
| No FDR correction on cell fraction tests | 04_cell_fractions.R | Added `p.adjust(method = "BH")` |
| Undefined `SeuratObj_ERpos_Young/Elderly` | 09_cellphonedb_prep.R | Properly subset using `AgeGroup` |
| Typo `CountData_YoungElderlyroWhole` | 09_cellphonedb_prep.R | Fixed variable name |
| Typo `_Macrophage` vs `_MonoMacro` | 03_cell_type_annotation.R | Check `CellTypeMinor` exists before use |
| Missing GSVA code | - | Created 06_run_gsva.R (pseudo-bulk) |
| Missing PROGENy code | - | Created 07_run_progeny.R (pseudo-bulk) |
| Missing WCSEA code | - | Created 08_run_wcsea.R |

### Rat Bulk RNA-seq (Zip)

| Issue | File | Fix |
|-------|------|-----|
| DESeq2 results not FDR-filtered | 03_deseq2.R | Filter by `padj < 0.05` not raw p-value |
| `input_dir1` vs `input_dir` variable | 01_alignment.sh | Changed to `INPUT_DIR` consistently |
| Missing SBATCH mail directives | All .sh files | Added `--mail-type=FAIL --mail-user=` |
| `annot.nkis` inappropriate for rat | 04_pam50_subtyping.R | Create minimal annotation data frame |
| HTSeq missing strand flag | 02_htseq_count.sh | Added `-s $STRAND` parameter |

## Reproducibility Gaps Filled

1. **scRNA GSVA**: Created pseudo-bulk approach aggregating per patient
2. **scRNA PROGENy**: Same pseudo-bulk approach
3. **WCSEA**: Implemented using indepthPathway package, using DE genes as input
```

**Step 2: Commit**

```bash
git add docs/plans/2026-02-17-integration-issues-log.md
git commit -m "docs: add integration issues log documenting all fixes"
```

---

## Final Verification

### Task 18: Run lint checks and verify structure

**Step 1: Verify directory structure**

```bash
ls -la analysis/04_human_scrnaseq/
ls -la analysis/05_rat_bulk_rnaseq/
```

Expected output shows all scripts and run files.

**Step 2: Check R scripts have proper headers**

```bash
head -5 analysis/04_human_scrnaseq/*.R
head -5 analysis/05_rat_bulk_rnaseq/*.R
```

**Step 3: Final commit**

```bash
git status
git log --oneline -10
```

---

**Plan complete and saved to `docs/plans/2026-02-17-code-integration-implementation.md`.**

**Two execution options:**

**1. Subagent-Driven (this session)** - I dispatch fresh subagent per task, review between tasks, fast iteration

**2. Parallel Session (separate)** - Open new session with executing-plans, batch execution with checkpoints

**Which approach?**