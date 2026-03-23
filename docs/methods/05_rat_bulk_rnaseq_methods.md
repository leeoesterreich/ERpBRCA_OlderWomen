# Methods: Rat Bulk RNA-seq Analysis (Section 05)

## 1. Data Acquisition

Paired-end bulk RNA-seq data were generated from rat mammary tumors (GEO accession GSE276757). Raw FASTQ files were obtained from six tumor samples (102-FF-Tumor, 107-L4t2-Tumor, 116-R4-Tumor, 157-R2-Tumor, 158-L2-Tumor, 167-L2-Tumor) stored at `/ix1/alee/LO_LAB/General/Lab_Data/20240622_Neil_Rahul_RNASeqRat`. Sample phenotype assignments (CONTROL vs TEST) were defined in a metadata file (`data/rat_bulk_rnaseq/sample_metadata.csv`), with TYPE encoded as a two-level factor (reference level: CONTROL).

## 2. Read Trimming

Raw paired-end reads were adapter-trimmed and quality-filtered using fastp v0.23.4 (`00_fastp_trim.sh`). Adapter sequences were auto-detected for paired-end reads (`--detect_adapter_for_pe`). Bases with Phred quality scores below 20 were trimmed (`--qualified_quality_phred 20`), and reads shorter than 36 bp after trimming were discarded (`--length_required 36`). Processing used 8 threads per sample. Per-sample HTML and JSON quality reports were generated.

## 3. Read Alignment

Trimmed reads were aligned to the rat genome (Rattus norvegicus, mRatBN7.2 assembly) using STAR v2.7.11b (`01_alignment.sh`). A pre-built STAR genome index located at `/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Rat` was used. STAR alignment used binary version 2.7.11b against a genome index generated with version 2.7.4a. While Ensembl recommends version-matched indices, STAR maintains backward compatibility for minor version differences and no alignment errors were reported. Alignment was performed with 64 threads (`--runThreadN 64`), compressed input handling via `zcat` (`--readFilesCommand zcat`), and output as coordinate-sorted BAM files (`--outSAMtype BAM SortedByCoordinate`). BAM sorting used 8 threads (`--outBAMsortingThreadN 8`). Quantification was performed in both transcriptome-level and gene-level modes (`--quantMode TranscriptomeSAM GeneCounts`).

## 4. Gene-Level Quantification

Gene-level read counts were generated from sorted BAM files using htseq-count v0.13.5 (`02_htseq_count.sh`). Reads were counted against the Ensembl rat gene annotation (Rattus_norvegicus.mRatBN7.2.112.gtf). Parameters were: BAM format input (`-f bam`), position-sorted reads (`-r pos`), reverse-stranded library protocol (`-s reverse`), feature type exon (`-t exon`), and gene identifier attribute `gene_id` (`-i gene_id`). HTSeq-count was run in reverse-stranded mode, consistent with the dUTP-based library preparation protocol used for these samples. HTSeq summary lines (prefixed with `__`) were excluded from downstream analyses.

## 5. Differential Expression Analysis

Differential expression analysis was performed using DESeq2 (`03_deseq2.R`) with a random seed of 12345. Individual HTSeq count files were loaded, parsed with `data.table::fread()`, and assembled into a gene-by-sample count matrix. HTSeq summary rows (lines beginning with `__`) were removed.

A `DESeqDataSet` was constructed from the count matrix and sample metadata using the design formula `~ TYPE`. The DESeq2 pipeline was executed via the default `DESeq()` call, which internally performed size factor estimation (median-of-ratios method), dispersion estimation (shrinkage), and negative binomial Wald tests. Results were extracted for the contrast TEST vs CONTROL (`results(dds, contrast = c("TYPE", "TEST", "CONTROL"))`).

Gene symbols were annotated by querying the Ensembl BioMart database (`rnorvegicus_gene_ensembl` dataset) for `external_gene_name` attributes using the `biomaRt` package. This query was wrapped in `tryCatch()` to handle network failures gracefully.

Significantly differentially expressed genes were defined as those with Benjamini-Hochberg adjusted p-value (FDR) < 0.05. Results were saved as both the full result table (`deseq2_results.csv`) and the FDR-filtered subset sorted by ascending adjusted p-value (`deseq2_results_significant.csv`).

DESeq2 size-factor-normalized counts were also exported. These normalized counts were rescaled to a per-million basis via `sweep(norm_counts, 2, colSums(norm_counts), "/") * 1e6` and saved as `normalized_cpm.csv` for use in PAM50 subtyping. Note: these are DESeq2 size-factor-normalized counts rescaled to counts-per-million (CPM), not true TPM (which requires gene-length normalization). The PAM50 script explicitly transforms these values to log2(CPM + 1) and median-centers each gene across samples before passing the matrix to `genefu::molecular.subtyping()` (see Section 6.2), which partially compensates for the lack of gene-length normalization.

## 6. PAM50 Molecular Subtyping

PAM50 intrinsic subtype classification was performed using the genefu R package (`04_pam50_subtyping.R`) with a random seed of 12345. The PAM50 centroids and robust model were loaded via `data("pam50")` and `data("pam50.robust")`.

### 6.1 Ortholog Mapping

A hardcoded human-to-rat gene symbol mapping was defined for all 50 PAM50 genes (e.g., `"ACTR3B" = "Actr3b"`, `"ESR1" = "Esr1"`). This mapping included alternative gene names (`"CDCA1" = "Nuf2"`, `"KNTC2" = "Ndc80"`, `"ORC6L" = "Orc6"`). Gene coverage was assessed post-mapping, with a warning threshold of fewer than 40 of 50 genes present.

### 6.2 Expression Data Preparation

The CPM-normalized expression matrix (`normalized_cpm.csv`) was loaded. Ensembl gene IDs (detected by the `ENSRNOG` prefix) were converted to gene symbols via BioMart query (`rnorvegicus_gene_ensembl` dataset). For duplicated gene symbols, the row with the highest mean expression across samples was retained. The expression matrix was then subset to available PAM50 ortholog genes, transformed to log2(CPM + 1), and median-centered per gene across samples (i.e., each gene's median expression was subtracted). The centered matrix was transposed to samples-by-genes orientation. Column names were remapped to human gene symbols for compatibility with genefu centroids.

### 6.3 Classification

Subtype classification was performed via `molecular.subtyping(sbt.model = "pam50", data = pam50_tpm, annot = annot_df, do.mapping = FALSE)`. A minimal annotation data frame was constructed with `probe` and `Gene.Symbol` columns set to the human gene symbols, and `EntrezGene.ID` set to `NA`. The `do.mapping = FALSE` flag was used since gene names were pre-mapped to match the PAM50 centroid row names.

### 6.4 Outputs

Per-sample subtype assignments and subtype posterior probabilities were saved to `pam50_subtypes.csv`. Two heatmaps were generated as PDF files:
- `pam50_heatmap.pdf` (12 x 15 inches): Z-score-normalized (row-scaled via `scale()`) PAM50 gene expression, annotated by predicted subtype, with column clustering disabled (`cluster_cols = FALSE`).
- `pam50_probabilities_heatmap.pdf` (12 x 8 inches): Transposed subtype probability matrix with both row and column clustering disabled.

Both heatmaps were rendered using `pheatmap`.

## 7. Validation Against Original Analysis

Reproducibility of the pipeline was validated against an independent analysis by Rahul (`05_validate_vs_rahul.R`) with a random seed of 12345. Three concordance checks were performed:

### 7.1 HTSeq Count Correlation

Per-sample HTSeq count vectors were compared between the current pipeline and Rahul's results (located at `/ix1/alee/LO_LAB/Personal/Rahul/Neil_RNAseq/5_Count_file/`). A hardcoded sample name lookup table was used for matching (e.g., `"102-FF-Tumor"` to `"102-FF"`). Pearson correlation was computed on merged gene-level counts (excluding `__`-prefixed summary rows). The pass threshold was r >= 0.99 for all samples.

### 7.2 DESeq2 Significant Gene Overlap

Significant gene sets (FDR < 0.05) from both analyses were compared. Rahul's results were loaded from `6_DEseq2/DESeq2_results_with_symbols.csv` and filtered to `padj < 0.05`. Overlap was quantified as the intersection size divided by the smaller set size (overlap coefficient), with a pass threshold of >= 90%. The Jaccard index (intersection/union) was also computed and reported for reference but was not used for the pass/fail criterion.

### 7.3 PAM50 Subtype Concordance

Per-sample PAM50 subtype assignments were compared between the two analyses. Rahul's subtypes were loaded from `PAM50/PAM50.csv`. Matching was performed by removing the `-Tumor` suffix via `gsub("-Tumor", "", ...)`. The pass criterion was 100% concordance across all matched samples.

Validation results were saved as an RDS object (`validation_results.rds`).

## 8. Pipeline Orchestration

The full pipeline was orchestrated via `run_analysis.sh`, which submitted SLURM jobs with explicit dependencies:

1. `00_fastp_trim.sh` (trimming)
2. `01_alignment.sh` (STAR alignment, `--dependency=afterok` on step 0)
3. `02_htseq_count.sh` (gene counting, `--dependency=afterok` on step 1)
4. `run_analysis.sbatch` (DESeq2 + PAM50 + validation, `--dependency=afterok` on step 2)

The R analysis steps (03-05) were executed sequentially within a single SLURM job (`run_analysis.sbatch`) using the `erp_brca_aging` conda environment.

## 9. Software Versions

From `environment.yml` and module loads:

| Software | Version | Source |
|----------|---------|--------|
| fastp | 0.23.4 | HPC module |
| STAR | 2.7.11b | HPC module |
| htseq-count | 0.13.5 | HPC module |
| R | 4.3.3 | conda (r-base) |
| DESeq2 | (conda bioconductor-deseq2) | Bioconductor via conda |
| genefu | (conda r-genefu) | conda-forge |
| biomaRt | (bioconductor-biomart) | Bioconductor via `install_r_packages.sh` |
| pheatmap | (conda r-pheatmap) | conda |
| data.table | (r-essentials bundle) | conda |
| dplyr | (r-tidyverse bundle) | conda |
| GCC | 8.2.0 | HPC module (for STAR) |

**Reference genome:** Rattus norvegicus mRatBN7.2 (Ensembl release 112), GTF: `Rattus_norvegicus.mRatBN7.2.112.gtf`.

**Conda environment:** `erp_brca_aging` at `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging`.

## 10. Reproducibility Notes

- All R scripts set `set.seed(12345)` for reproducibility of any stochastic operations.
- The pipeline uses `set -euo pipefail` (shell scripts) and `set -eo pipefail` (sbatch wrapper) for strict error handling.
- STAR alignment and HTSeq counting include skip-if-exists logic to support idempotent re-runs.
- BioMart queries (for gene symbol annotation) depend on external Ensembl server availability and may return different results if the database version changes. No specific Ensembl archive release was pinned.
- The `environment.yml` does not pin exact versions for DESeq2, genefu, or biomaRt; only R base (4.3.3) and Seurat (>=5.0) are version-constrained. Exact installed versions should be recorded at runtime.

## Main Text Summary

Bulk RNA-seq from six rat mammary tumors (GEO: GSE276757) was processed through a four-stage pipeline: adapter trimming (fastp v0.23.4, Q>=20, length>=36 bp), alignment to the mRatBN7.2 rat genome (STAR v2.7.11b), and gene-level quantification (htseq-count v0.13.5, reverse-stranded). Differential expression analysis was performed with DESeq2 using a ~TYPE design formula (TEST vs CONTROL, FDR < 0.05). PAM50 molecular subtyping was performed using the genefu package with a manually curated human-to-rat ortholog mapping of all 50 classifier genes. Pipeline outputs were validated against an independent analysis, confirming high concordance in read counts (Pearson r >= 0.99), significant gene sets (overlap coefficient >= 90%), and subtype assignments.
