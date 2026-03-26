# Data Directory

This directory contains input data and download scripts for the analyses.

## GEO Accessions

| Dataset | GEO Accession | Description |
|---------|---------------|-------------|
| Human Bulk RNA-seq | [GSE276755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276755) | 168 samples (83 tumor, 85 adjacent normal) |
| Rat snRNA-seq | [GSE276758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276758) | 6 samples (3 Young, 3 Aged) |
| Rat WES | [GSE276759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276759) | Whole exome sequencing of rat tumors |
| Human scRNA-seq | Xu et al. 2024 atlas | HR+ treatment-naive breast tumor atlas (~115K cells) |

## Download Instructions

Run the master download script:

```bash
bash download_all.sh
```

Or download individual datasets:

```bash
bash human_bulk_rnaseq/download.sh
bash rat_snrnaseq/download.sh
bash rat_wes/download.sh
```

The human scRNA-seq reanalysis uses the Xu et al. 2024 primary breast tumor atlas. See `human_scrnaseq/README.md` for details on obtaining the atlas data.

## Directory Structure

```
data/
├── download_all.sh               # Master download script
├── gmt/                          # Custom GMT gene sets for pathway analyses
├── human_bulk_rnaseq/
│   ├── download.sh
│   └── external/                 # Gene annotation files
├── human_scrnaseq/
│   ├── README.md                 # Xu et al. atlas data provenance
│   ├── metadata/
│   └── raw/                      # Downloaded count matrices
├── rat_bulk_rnaseq/
│   ├── README.md
│   └── sample_metadata.csv
├── rat_snrnaseq/
│   ├── download.sh
│   ├── metadata/
│   └── raw/                      # 10X filtered matrices per sample
└── rat_wes/
    └── download.sh
```

## MICA Input Data

The MICA concordance analysis (`analysis/01_human_bulk_rnaseq/mica/`) uses pre-computed pathway scores and expression matrices derived from TCGA, METABRIC, and SCAN-B. These datasets cannot be redistributed due to their respective data access agreements and must be obtained from their original sources (see main README for links).
