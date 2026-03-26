# Chronic Inflammation and Hormone Disposition in Breast Cancer of Older Women

Companion analysis repository for the Nature Aging study on aging-associated estrogen signaling and inflammatory remodeling in ER-positive breast cancer.

## Overview

This repository contains analysis code associated with the study. The codebase is organized by assay or analysis module, with separate workflows for human bulk RNA-seq, rat single-nucleus RNA-seq, rat whole-exome sequencing, human single-cell RNA-seq reanalysis, rat bulk RNA-seq, spatial biopsy analyses, and organoid single-cell RNA-seq.

## Key Repository Layout

```text
.
├── analysis/
│   ├── 01_human_bulk_rnaseq/      # Human bulk RNA-seq preprocessing, GSVA, PROGENy, correlations, figures
│   ├── 02_rat_snrnaseq/           # Rat snRNA-seq QC, integration, annotation, DE, DA
│   ├── 03_rat_wes/                # Rat WES parsing, COSMIC signatures, oncoplot generation
│   ├── 04_human_scrnaseq/         # Human scRNA-seq reanalysis and macrophage-focused downstream analyses
│   ├── 05_rat_bulk_rnaseq/        # Rat bulk RNA-seq alignment, counting, DESeq2, PAM50 subtyping
│   ├── 06_spatial_biopsies/       # Spatial biopsy cytokine and pathway analyses
│   └── 07_organoid_single_cell/   # Organoid single-cell RNA-seq analysis
├── data/
│   ├── gmt/                       # Custom GMT gene sets used by pathway analyses
│   ├── human_bulk_rnaseq/
│   ├── human_scrnaseq/
│   ├── rat_bulk_rnaseq/
│   ├── rat_snrnaseq/
│   └── rat_wes/
├── environment.yml                # Conda environment specification
└── sbatch_run_all.sh              # Cluster launcher for local HPC use
```

## Analysis Modules

### 01. Human bulk RNA-seq

`analysis/01_human_bulk_rnaseq/` runs a five-step workflow: preprocessing of the human bulk expression matrix and sample annotation, GSVA for estrogen-related gene sets, PROGENy pathway activity inference, Spearman correlation analysis between selected genes and pathway scores, and export of correlation bubble plots.

Outputs are written to `analysis/01_human_bulk_rnaseq/outputs/`.

### 02. Rat snRNA-seq

`analysis/02_rat_snrnaseq/` contains scripts for QC and filtering of six nuclei datasets, SCTransform normalization, Harmony integration, clustering, marker-based annotation, pseudobulk differential expression by cell type, and differential abundance testing with `speckle::propeller`.

Key outputs are written to `analysis/02_rat_snrnaseq/outputs/`.

### 03. Rat WES

`analysis/03_rat_wes/` parses VEP output, summarizes COSMIC mutational signatures, and generates oncoplot-ready outputs and figures. The executable cluster entry point in this directory is `run_analysis.sbatch`.

### 04. Human scRNA-seq reanalysis

`analysis/04_human_scrnaseq/` includes a GEO download helper (`00a_download_geo.sh`), a preprocessing script for an HR-positive, treatment-naive subset of the Xu et al. breast tumor atlas (`00b_preprocess_seurat.R`), and downstream annotation, cell-fraction analysis, gene-expression visualization, GSVA, PROGENy, WCSEA, CellPhoneDB preparation, CellChat analysis, and macrophage-focused differential/pathway analyses.

### 05. Rat bulk RNA-seq

`analysis/05_rat_bulk_rnaseq/` contains trimming, alignment, HTSeq counting, DESeq2 differential expression, and PAM50 molecular subtyping scripts. The shell launcher submits SLURM jobs for the preprocessing steps and then runs the downstream R analysis.

### 06. Spatial biopsies

`analysis/06_spatial_biopsies/` generates spatial plots for cytokines including IL4, IL10, IL13, TGFB1, TGFB2, and TGFB3, along with pathway-enrichment and preranked GSEA summaries focused on IL/TGF-related immune programs. Input data must be obtained separately (see Data Availability).

### 07. Organoid single-cell RNA-seq

`analysis/07_organoid_single_cell/` contains the organoid workflow from Cell Ranger outputs through QC, preprocessing, pseudobulk differential expression, pathway analysis, cell-cycle scoring, heterogeneity analysis, and proliferation analysis. This module also contains its own `README.md` with assay-specific details.

## Data Availability

The primary study dataset is available at GEO under accession [GSE276755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276755).

Additional datasets referenced in the repository include [GSE276758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276758) for rat snRNA-seq, [GSE276759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276759) for rat whole-exome sequencing, and the Xu et al. breast tumor atlas for human scRNA-seq reanalysis.

The MICA concordance analysis (`analysis/01_human_bulk_rnaseq/mica/`) uses pre-computed pathway scores and expression matrices derived from TCGA, METABRIC, and SCAN-B. These datasets cannot be redistributed due to their respective data access agreements and must be obtained directly from their original sources:
- **TCGA**: [GDC Data Portal](https://portal.gdc.cancer.gov/) (dbGaP access required)
- **METABRIC**: [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric) or [EGA](https://ega-archive.org/studies/EGAS00000000083)
- **SCAN-B**: [GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) (expression) and [GSE81540](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81540) (clinical)

The helper script `data/download_all.sh` downloads the GEO-hosted study datasets (`GSE276755`, `GSE276758`, `GSE276759`). Some workflows also reference local institutional paths for source data that are not mirrored in this repository.

## Environment and Reproducibility

A base conda environment is defined in `environment.yml`:

```bash
conda env create -f environment.yml
conda activate erp_brca_aging
```

The repository also includes SLURM-oriented launchers such as module-level `run_analysis.sbatch` scripts and `sbatch_run_all.sh`.

In practice, the module-level launchers are the clearest entry points for rerunning analyses. Scripts assume a SLURM environment. Some data paths (reference genomes, raw FASTQs) are configured via environment variables and will need to be set for your local compute environment.

Example entry points:

```bash
bash data/download_all.sh
bash analysis/01_human_bulk_rnaseq/run_analysis.sh
bash analysis/02_rat_snrnaseq/run_analysis.sh
sbatch analysis/03_rat_wes/run_analysis.sbatch
bash analysis/05_rat_bulk_rnaseq/run_analysis.sh
bash analysis/06_spatial_biopsies/run_analysis.sh
sbatch analysis/07_organoid_single_cell/run_analysis.sbatch
```

For human scRNA-seq, `analysis/04_human_scrnaseq/run_analysis.sh` submits the corresponding SLURM job.

## Issues

If you encounter problems running these analyses or have questions about the code, please [open a GitHub issue](https://github.com/leeoesterreich/ERpBRCA_OlderWomen/issues) and someone from the lab will respond.

## Citation

Please cite the accompanying Nature Aging article and this repository when using these materials.

## License

This repository is distributed under the terms of the `LICENSE` file.
