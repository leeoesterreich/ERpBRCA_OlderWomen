# Data Directory

This directory contains all input data for the analyses.

## GEO Accessions

| Dataset | GEO Accession | Description |
|---------|---------------|-------------|
| Human Bulk RNA-seq | GSE276755 | 168 samples (83 tumor, 85 adjacent normal) |
| Rat snRNA-seq | GSE276758 | 6 samples (3 Young, 3 Aged) |
| Rat WES | GSE276759 | Whole exome sequencing of rat tumors |

## Download Instructions

Run the master download script:

```bash
./download_all.sh
```

Or download individual datasets:

```bash
./human_bulk_rnaseq/download.sh
./rat_snrnaseq/download.sh
./rat_wes/download.sh
```

## Directory Structure

```
data/
├── human_bulk_rnaseq/
│   ├── raw/                    # Count and TPM matrices
│   └── external/               # Gene annotation files
├── rat_snrnaseq/
│   ├── raw/                    # 10X filtered matrices per sample
│   └── metadata/               # Sample annotations
└── rat_wes/
    └── raw/                    # VEP and CNVKit outputs
```
