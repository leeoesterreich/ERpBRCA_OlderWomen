# Rat Bulk RNA-seq Pipeline Test Design

**Date:** 2026-02-23
**Status:** Approved
**Goal:** Full re-run of rat bulk RNA-seq pipeline with validation against Rahul's original results

## Context

The code integration (Feb 17) created refactored scripts in `analysis/05_rat_bulk_rnaseq/`. Raw data is now available locally. This test will:
1. Run the full pipeline from raw FASTQs
2. Validate outputs match Rahul's original analysis

## Data

### Samples
| Sample | Group | Description |
|--------|-------|-------------|
| 102-FF-Tumor | TEST | Young rat tumor |
| 107-L4t2-Tumor | TEST | Young rat tumor |
| 116-R4-Tumor | TEST | Young rat tumor |
| 157-R2-Tumor | CONTROL | Elderly rat tumor |
| 158-L2-Tumor | CONTROL | Elderly rat tumor |
| 167-L2-Tumor | CONTROL | Elderly rat tumor |

### Paths
- **Raw FASTQs:** `/ix1/alee/LO_LAB/General/Lab_Data/20240622_Neil_Rahul_RNASeqRat/`
- **Reference genome:** `/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Rat/`
- **Rahul's outputs:** `/ix1/alee/LO_LAB/Personal/Rahul/Neil_RNAseq/`
- **Metadata:** `data/rat_bulk_rnaseq/sample_metadata.csv`

## Pipeline

```
Raw FASTQs (.fastq.gz)
    ↓
00_fastp_trim.sh (NEW)
    ↓
01_alignment.sh (STAR → sorted BAM)
    ↓
02_htseq_count.sh (gene counts)
    ↓
03_deseq2.R (differential expression)
    ↓
04_pam50_subtyping.R (molecular subtyping)
    ↓
05_validate_vs_rahul.R (NEW - comparison)
```

## Scripts to Create/Modify

### NEW: 00_fastp_trim.sh
- Input: Raw compressed FASTQs
- Output: Trimmed FASTQs + QC reports
- Tool: fastp (auto-detect adapters)

### MODIFY: 01_alignment.sh
- Update GENOME_DIR to Rahul's reference
- Update INPUT_DIR to trimmed output
- Add `--readFilesCommand zcat` for compressed input

### MODIFY: 02_htseq_count.sh
- Update paths
- Confirm strand setting (reverse for typical Illumina)

### MODIFY: 03_deseq2.R
- Update paths to our count files
- Use sample_metadata.csv for sample info

### MODIFY: 04_pam50_subtyping.R
- Update paths
- Generate normalized TPM for input

### NEW: 05_validate_vs_rahul.R
- Compare HTSeq counts (Pearson correlation)
- Compare DESeq2 results (significant gene overlap)
- Compare PAM50 calls (exact match)
- Generate validation_report.txt

## Validation Criteria

| Check | Pass Criteria |
|-------|---------------|
| HTSeq counts | Pearson r ≥ 0.99 per sample |
| DESeq2 significant genes | ≥90% overlap at FDR<0.05 |
| PAM50 subtypes | 6/6 identical calls |

## Output Structure

```
analysis/05_rat_bulk_rnaseq/outputs/
├── trimmed/
│   ├── 102-FF-Tumor_R1_trimmed.fastq.gz
│   ├── 102-FF-Tumor_R2_trimmed.fastq.gz
│   ├── 102-FF-Tumor_fastp.html
│   └── ...
├── aligned/
│   ├── 102-FF-Tumor_Aligned.sortedByCoord.out.bam
│   └── ...
├── counts/
│   ├── 102-FF-Tumor.txt
│   └── ...
├── deseq2_results.csv
├── deseq2_results_significant.csv
├── pam50_subtypes.csv
├── pam50_heatmap.pdf
└── validation_report.txt
```

## Runtime Estimate

| Step | Time |
|------|------|
| Trimming (6 samples) | ~30 min |
| Alignment (6 samples) | ~4-5 hours |
| HTSeq counting | ~1 hour |
| DESeq2 + PAM50 | ~10 min |
| **Total** | **~6-7 hours** |

## SLURM Resources

- **Alignment:** 64 cores, 1 day limit
- **HTSeq:** 8 cores, 12 hours
- **R analysis:** 4 cores, 32GB RAM, 4 hours
