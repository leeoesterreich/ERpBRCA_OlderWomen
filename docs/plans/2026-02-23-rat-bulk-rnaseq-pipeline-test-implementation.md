# Rat Bulk RNA-seq Pipeline Test Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Run the full rat bulk RNA-seq pipeline from raw FASTQs and validate outputs match Rahul's original analysis.

**Architecture:** Sequential pipeline: fastp trimming → STAR alignment → HTSeq counting → DESeq2 DE analysis → PAM50 subtyping → validation comparison. Each step produces outputs consumed by the next.

**Tech Stack:** Bash (fastp, STAR, HTSeq), R (DESeq2, genefu, biomaRt), SLURM job scheduler

---

## Task 1: Create fastp trimming script

**Files:**
- Create: `analysis/05_rat_bulk_rnaseq/00_fastp_trim.sh`

**Step 1: Create the trimming script**

```bash
#!/bin/bash
#SBATCH --job-name=rat_fastp
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 4:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/fastp_%j.out
#SBATCH --error=logs/fastp_%j.err

# Rat bulk RNA-seq trimming with fastp
# Inputs: Raw FASTQ files
# Outputs: Trimmed FASTQ files + QC reports

set -euo pipefail

module purge
module load fastp/0.23.4

# Paths
RAW_DIR="/ix1/alee/LO_LAB/General/Lab_Data/20240622_Neil_Rahul_RNASeqRat"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/outputs/trimmed"

mkdir -p "$OUTPUT_DIR"
mkdir -p "${SCRIPT_DIR}/logs"

echo "=== fastp Trimming ==="
echo "Input: $RAW_DIR"
echo "Output: $OUTPUT_DIR"

# Process each sample
for read1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
    [ -f "$read1" ] || continue

    # Extract sample name (e.g., 102-FF-Tumor from 102-FF-Tumor_R1_001.fastq.gz)
    sample=$(basename "$read1" _R1_001.fastq.gz)
    read2="${RAW_DIR}/${sample}_R2_001.fastq.gz"

    if [ ! -f "$read2" ]; then
        echo "Warning: R2 not found for $sample, skipping"
        continue
    fi

    echo "Processing: $sample"

    fastp \
        -i "$read1" \
        -I "$read2" \
        -o "${OUTPUT_DIR}/${sample}_R1_trimmed.fastq.gz" \
        -O "${OUTPUT_DIR}/${sample}_R2_trimmed.fastq.gz" \
        --html "${OUTPUT_DIR}/${sample}_fastp.html" \
        --json "${OUTPUT_DIR}/${sample}_fastp.json" \
        --thread 8 \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 36

    echo "Completed: $sample"
done

echo "=== Trimming complete ==="
echo "Output files:"
ls -la "$OUTPUT_DIR"
```

**Step 2: Make executable and verify syntax**

Run: `chmod +x analysis/05_rat_bulk_rnaseq/00_fastp_trim.sh && bash -n analysis/05_rat_bulk_rnaseq/00_fastp_trim.sh`
Expected: No output (syntax OK)

**Step 3: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/00_fastp_trim.sh
git commit -m "feat(05_rat_bulk): add fastp trimming script"
```

---

## Task 2: Update alignment script with correct paths

**Files:**
- Modify: `analysis/05_rat_bulk_rnaseq/01_alignment.sh`

**Step 1: Read current script**

Read the existing script to understand current structure.

**Step 2: Update the script with correct paths**

```bash
#!/bin/bash
#SBATCH --job-name=rat_alignment
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00
#SBATCH --mem=128G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/alignment_%j.out
#SBATCH --error=logs/alignment_%j.err

# Rat bulk RNA-seq alignment with STAR
# Inputs: Trimmed FASTQ files from fastp
# Outputs: Sorted BAM files

set -euo pipefail

module purge
module load gcc/8.2.0
module load star/2.7.11b

# Paths - UPDATED
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/outputs/trimmed"
OUTPUT_DIR="${SCRIPT_DIR}/outputs/aligned"
GENOME_DIR="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Rat"
THREADS=64

mkdir -p "$OUTPUT_DIR"
mkdir -p "${SCRIPT_DIR}/logs"

echo "=== STAR Alignment ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Genome: $GENOME_DIR"

# Verify genome index exists
if [ ! -f "${GENOME_DIR}/SAindex" ]; then
    echo "ERROR: STAR index not found at $GENOME_DIR"
    exit 1
fi

# Align each sample
for read1 in "$INPUT_DIR"/*_R1_trimmed.fastq.gz; do
    [ -f "$read1" ] || continue

    sample=$(basename "$read1" _R1_trimmed.fastq.gz)
    read2="${INPUT_DIR}/${sample}_R2_trimmed.fastq.gz"
    output_prefix="${OUTPUT_DIR}/${sample}_"

    if [ ! -f "$read2" ]; then
        echo "Warning: R2 not found for $sample, skipping"
        continue
    fi

    # Skip if already aligned
    if [ -f "${output_prefix}Aligned.sortedByCoord.out.bam" ]; then
        echo "Skipping $sample (already aligned)"
        continue
    fi

    echo "Aligning: $sample"
    STAR --genomeDir "$GENOME_DIR" \
         --readFilesIn "$read1" "$read2" \
         --readFilesCommand zcat \
         --runThreadN $THREADS \
         --outFileNamePrefix "$output_prefix" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts \
         --outBAMsortingThreadN 8

    echo "Completed: $sample"
done

echo "=== Alignment complete ==="
echo "Output files:"
ls -la "$OUTPUT_DIR"/*.bam 2>/dev/null || echo "No BAM files yet"
```

**Step 3: Verify syntax**

Run: `bash -n analysis/05_rat_bulk_rnaseq/01_alignment.sh`
Expected: No output (syntax OK)

**Step 4: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/01_alignment.sh
git commit -m "fix(05_rat_bulk): update alignment paths and add zcat for compressed input"
```

---

## Task 3: Update HTSeq counting script

**Files:**
- Modify: `analysis/05_rat_bulk_rnaseq/02_htseq_count.sh`

**Step 1: Update the script**

```bash
#!/bin/bash
#SBATCH --job-name=rat_htseq
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/htseq_%j.out
#SBATCH --error=logs/htseq_%j.err

# HTSeq read counting for rat bulk RNA-seq
# Inputs: Aligned BAM files
# Outputs: Gene count files

set -euo pipefail

module purge
module load htseq/0.13.5

# Paths - UPDATED
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/outputs/aligned"
OUTPUT_DIR="${SCRIPT_DIR}/outputs/counts"
GTF_FILE="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Rat/Rattus_norvegicus.mRatBN7.2.112.gtf"

# Strand setting - reverse for typical Illumina TruSeq
STRAND="reverse"

mkdir -p "$OUTPUT_DIR"

echo "=== HTSeq Count ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "GTF: $GTF_FILE"
echo "Strand: $STRAND"

for bam_file in "$INPUT_DIR"/*Aligned.sortedByCoord.out.bam; do
    [ -f "$bam_file" ] || continue

    sample=$(basename "$bam_file" _Aligned.sortedByCoord.out.bam)
    output_file="${OUTPUT_DIR}/${sample}.txt"

    # Skip if already counted
    if [ -f "$output_file" ]; then
        echo "Skipping $sample (already counted)"
        continue
    fi

    echo "Counting: $sample"
    htseq-count \
        -f bam \
        -r pos \
        -s "$STRAND" \
        -t exon \
        -i gene_id \
        "$bam_file" \
        "$GTF_FILE" > "$output_file"

    echo "Created: $output_file"
done

echo "=== HTSeq complete ==="
echo "Output files:"
ls -la "$OUTPUT_DIR"
```

**Step 2: Verify syntax**

Run: `bash -n analysis/05_rat_bulk_rnaseq/02_htseq_count.sh`
Expected: No output (syntax OK)

**Step 3: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/02_htseq_count.sh
git commit -m "fix(05_rat_bulk): update HTSeq paths and strand setting"
```

---

## Task 4: Update DESeq2 script

**Files:**
- Modify: `analysis/05_rat_bulk_rnaseq/03_deseq2.R`

**Step 1: Read current script to understand structure**

**Step 2: Update the script with correct paths**

Key changes:
- Update paths to use our outputs
- Use sample_metadata.csv from data directory
- Ensure FDR filtering is applied

```r
#!/usr/bin/env Rscript
# analysis/05_rat_bulk_rnaseq/03_deseq2.R
# Differential expression analysis for rat bulk RNA-seq

set.seed(12345)

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(biomaRt)
  library(dplyr)
})

# Define paths
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir)) script_dir <- getwd()

project_root <- normalizePath(file.path(script_dir, "../.."))
output_dir <- file.path(script_dir, "outputs")
count_dir <- file.path(output_dir, "counts")
data_dir <- file.path(project_root, "data/rat_bulk_rnaseq")

cat("=== DESeq2 Analysis ===\n")
cat("Count dir:", count_dir, "\n")
cat("Data dir:", data_dir, "\n")

# Step 1: Load count data
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

# Step 2: Load sample metadata
cat("\nStep 2: Loading sample metadata...\n")

metadata_file <- file.path(data_dir, "sample_metadata.csv")
coldata <- read.csv(metadata_file)
rownames(coldata) <- coldata$SAMPLE
coldata$TYPE <- factor(coldata$TYPE, levels = c("CONTROL", "TEST"))

# Match sample order
common_samples <- intersect(colnames(countdata), rownames(coldata))
countdata <- countdata[, common_samples]
coldata <- coldata[common_samples, ]

cat("  Samples:", nrow(coldata), "\n")
cat("  Groups:", paste(table(coldata$TYPE), collapse = " vs "), "\n")

# Step 3: Run DESeq2
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

# Step 4: Add gene symbols
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
  cat("  Warning: Could not connect to biomaRt\n")
  res_df$external_gene_name <- NA
})

# Step 5: Save results
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
```

**Step 3: Verify R syntax**

Run: `Rscript -e "parse('analysis/05_rat_bulk_rnaseq/03_deseq2.R')"`
Expected: No errors

**Step 4: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/03_deseq2.R
git commit -m "fix(05_rat_bulk): update DESeq2 paths and use sample_metadata.csv"
```

---

## Task 5: Update PAM50 subtyping script

**Files:**
- Modify: `analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R`

**Step 1: Update with correct paths and input handling**

Key changes:
- Read normalized TPM from DESeq2 output
- Update output paths

**Step 2: Verify R syntax**

Run: `Rscript -e "parse('analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R')"`

**Step 3: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R
git commit -m "fix(05_rat_bulk): update PAM50 paths to use DESeq2 output"
```

---

## Task 6: Create validation script

**Files:**
- Create: `analysis/05_rat_bulk_rnaseq/05_validate_vs_rahul.R`

**Step 1: Create the validation script**

```r
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
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir)) script_dir <- getwd()
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
```

**Step 2: Verify R syntax**

Run: `Rscript -e "parse('analysis/05_rat_bulk_rnaseq/05_validate_vs_rahul.R')"`

**Step 3: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/05_validate_vs_rahul.R
git commit -m "feat(05_rat_bulk): add validation script to compare with Rahul's results"
```

---

## Task 7: Update run scripts

**Files:**
- Modify: `analysis/05_rat_bulk_rnaseq/run_analysis.sh`
- Modify: `analysis/05_rat_bulk_rnaseq/run_analysis.sbatch`

**Step 1: Update run_analysis.sh to chain jobs**

```bash
#!/bin/bash
# Run rat bulk RNA-seq analysis pipeline
# Usage: bash run_analysis.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p logs outputs/trimmed outputs/aligned outputs/counts

echo "=== Rat Bulk RNA-seq Pipeline ==="
echo "Directory: $SCRIPT_DIR"
echo ""

# Step 0: Trimming
echo "Step 0: Submitting trimming job..."
JOB0=$(sbatch --parsable 00_fastp_trim.sh)
echo "  Job ID: $JOB0"

# Step 1: Alignment (depends on trimming)
echo "Step 1: Submitting alignment job..."
JOB1=$(sbatch --parsable --dependency=afterok:$JOB0 01_alignment.sh)
echo "  Job ID: $JOB1"

# Step 2: HTSeq counting (depends on alignment)
echo "Step 2: Submitting HTSeq job..."
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 02_htseq_count.sh)
echo "  Job ID: $JOB2"

# Step 3-5: R analysis (depends on counting)
echo "Step 3-5: Submitting R analysis job..."
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 run_analysis.sbatch)
echo "  Job ID: $JOB3"

echo ""
echo "=== Jobs submitted ==="
echo "Monitor with: squeue -u \$USER"
echo "Chain: $JOB0 -> $JOB1 -> $JOB2 -> $JOB3"
```

**Step 2: Update run_analysis.sbatch**

```bash
#!/bin/bash
#SBATCH --job-name=rat_bulk_R
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/R_analysis_%j.out
#SBATCH --error=logs/R_analysis_%j.err

# Rat bulk RNA-seq R analysis
set -euo pipefail

module purge
module load r/4.2.0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== DESeq2 Analysis ==="
Rscript 03_deseq2.R

echo "=== PAM50 Subtyping ==="
Rscript 04_pam50_subtyping.R

echo "=== Validation ==="
Rscript 05_validate_vs_rahul.R

echo "=== Analysis complete ==="
```

**Step 3: Commit**

```bash
git add analysis/05_rat_bulk_rnaseq/run_analysis.sh analysis/05_rat_bulk_rnaseq/run_analysis.sbatch
git commit -m "fix(05_rat_bulk): update run scripts with job chaining"
```

---

## Task 8: Submit pipeline and monitor

**Step 1: Submit the pipeline**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/05_rat_bulk_rnaseq && bash run_analysis.sh`

**Step 2: Monitor jobs**

Run: `squeue -u $USER`

**Step 3: Check logs as jobs complete**

Run: `tail -f logs/*.out`

---

## Task 9: Review validation results

**Step 1: Check validation report**

After pipeline completes, run:
```bash
cat analysis/05_rat_bulk_rnaseq/outputs/validation_results.rds
# Or in R:
Rscript -e "results <- readRDS('analysis/05_rat_bulk_rnaseq/outputs/validation_results.rds'); print(results)"
```

**Step 2: If all checks pass, commit final state**

```bash
git add -A
git commit -m "feat(05_rat_bulk): complete pipeline test - all validations passed"
```

**Step 3: If checks fail, investigate and iterate**

Review logs and outputs to identify discrepancies.
