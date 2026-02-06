# Chronic Inflammation and Hormone Disposition in Breast Cancer of Older Women

[![DOI](https://img.shields.io/badge/DOI-10.1101%2FXXXX-blue)](https://doi.org/10.1101/XXXX)
[![GEO](https://img.shields.io/badge/GEO-GSE276755-green)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276755)

**Carleton et al. 2024**

> Chronic inflammation and hormone disposition promote a tumor-permissive locale for breast cancer in older women

## Overview

This repository contains all analysis code for the study examining how the aged breast tumor microenvironment promotes ER+ breast cancer development through chronic inflammation and estrogen signaling.

### Key Findings

- Aged F344 rats develop tumors faster than younger counterparts
- Tumors upregulate HSD17B7 to convert E1 to E2 locally
- Tumor-associated macrophages integrate E2 and chemokine signaling
- Age-related immune dysfunction associated with chronic inflammation

## Repository Structure

```
├── analysis/
│   ├── 01_human_bulk_rnaseq/    # Human ER+ breast cancer RNA-seq
│   ├── 02_rat_snrnaseq/         # Rat mammary tumor snRNA-seq
│   └── 03_rat_wes/              # Rat whole exome sequencing
├── data/                         # Data download scripts (GEO)
├── figures/                      # Generated figures
├── results/                      # Analysis outputs
├── scripts/                      # Pipeline scripts
└── docs/                         # Documentation
```

## Quick Start

### 1. Setup Environment

```bash
# Clone repository
git clone https://github.com/leeoesterreich/ERpBRCA_OlderWomen.git
cd ERpBRCA_OlderWomen

# Create conda environment
conda env create -f environment.yml
conda activate erp_brca_aging
```

### 2. Download Data

```bash
# Download all data from GEO
bash data/download_all.sh
```

### 3. Run Analysis

```bash
# Run complete pipeline
bash scripts/run_all.sh

# Or run individual analyses
bash analysis/01_human_bulk_rnaseq/run_analysis.sh
bash analysis/02_rat_snrnaseq/run_analysis.sh
bash analysis/03_rat_wes/run_analysis.sh
```

## Data Availability

| Dataset | GEO Accession | Description |
|---------|---------------|-------------|
| Human Bulk RNA-seq | [GSE276755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276755) | 168 samples |
| Rat snRNA-seq | [GSE276758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276758) | 6 samples |
| Rat Bulk RNA-seq | [GSE276757](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276757) | - |
| Rat WES | [GSE276759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276759) | - |

## Methods Summary

### Human Bulk RNA-seq
- GSVA for estrogen pathway scoring
- PROGENy for pathway activity inference
- Spearman correlations with FDR correction

### Rat snRNA-seq
- Seurat v5 + Harmony integration
- DoubletFinder for doublet removal
- Differential expression (Wilcoxon + BH-FDR)
- Differential abundance (propeller)

## Translational Relevance

Paradoxically, ER+ breast cancer is an age-related disease, with the peak incidence occurring around the age of 67-70 years old, despite low circulating levels of estrogens. Here, we show that the unique features of the aged breast TME permit a conducive locale, including high local E2 and an immune dysfunctional environment due to chronic inflammation, for tumor development and growth. Our findings suggest that tumor-associated macrophages integrate the estrogen and chemokine-driven inflammation, polarize toward an immunosuppressive, and drive T cell exclusion. Our results provide evidence to support pharmacologically targeting immune components to re-polarize the aged tumor-immune microenvironment.

## Citation

```bibtex
@article{carleton2024chronic,
  title={Chronic inflammation and hormone disposition promote a tumor-permissive locale for breast cancer in older women},
  author={Carleton, Neil and others},
  journal={bioRxiv},
  year={2024},
  doi={10.1101/XXXX}
}
```

## License

This project is licensed under the terms of the [LICENSE](LICENSE) file.

## Contact

- Neil Carleton - [email]
- Lee-Oesterreich Lab - University of Pittsburgh
