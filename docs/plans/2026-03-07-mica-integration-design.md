# MICA Integration Design

**Date:** 2026-03-07
**Author:** Claude Code
**Status:** Approved

## Overview

Integrate Jian Zou's MICA (Mutual Information Concordance Analysis) code into the ERpBRCA_OlderWomen unified analysis pipeline. This enables reproduction of manuscript Figures 4C, 4H, 5C, and 5E.

## Background

MICA detects biomarkers with concordant multi-class expression patterns across multiple omics studies (TCGA, METABRIC, SCAN-B) using mutual information. The manuscript uses MICA for:

1. **Estrogen pathway analysis** (Fig 4) - EstroGene signatures, HSD17B enzymes
2. **Inflammation/immune analysis** (Fig 5) - Inflammatory pathways, immune cell deconvolution

## Source Code

Jian's original scripts in `Jian_MICA/`:
- `04_EstroGene_MICA[1].R` - Full data prep + GSVA + MICA pipeline
- `04_EstroGene_MICA.R` - Minimal MICA execution
- `03_240408_figure_update.R` - Figure generation (4C, 4H, 5C, 5E)

## Directory Structure

```
ERpBRCA_OlderWomen/
└── analysis/
    └── 01_human_bulk_rnaseq/
        └── mica/
            ├── 01_prep_data.R         # Load + filter TCGA/METABRIC/SCAN-B
            ├── 02_run_gsva.R          # GSVA on EstroGene + MSigDB pathways
            ├── 03_run_mica.R          # mica.full() permutation analysis
            ├── 04_fig4_estrogen.R     # Figs 4C (heatmap), 4H (boxplots)
            ├── 05_fig5_immune.R       # Figs 5C (heatmap), 5E (boxplots)
            ├── run_all.sbatch         # SLURM batch script
            ├── data/
            │   └── input/             # Symlinks to source data
            ├── outputs/               # .RData intermediate files
            ├── figures/               # PDF + PNG outputs
            └── logs/                  # SLURM logs
```

## Data Flow

```
data/input/ → 01_prep_data.R → outputs/01_*.RData
              02_run_gsva.R  → outputs/02_*.RData
              03_run_mica.R  → outputs/03_mica_results.RData
              04_fig4_*.R    → figures/fig4c_*, fig4h_*
              05_fig5_*.R    → figures/fig5c_*, fig5e_*
```

## Script Specifications

### 01_prep_data.R

**Purpose:** Load and filter expression data for ER+/HER2- samples

**Inputs:**
- `ilc_clean_data.RData` (TCGA, METABRIC, SCAN-B objects)
- `METABRIC_Clinical_Info.csv`
- `TCGA_Key.csv`

**Logic:**
1. Load BCdata objects (scanb, metabric, tcga)
2. Filter to ER+/HER2- samples
3. Create age groups:
   - Standard: Young (35-45), Middle-Aged (55-69), Elderly (70+)
   - PostM: Early (55-59), Middle (60-69), Elderly (70+)
4. Normalize expression (TMM for METABRIC microarray)

**Outputs:**
- `outputs/01_tcga_filtered.RData`
- `outputs/01_metabric_filtered.RData`
- `outputs/01_scanb_filtered.RData`

### 02_run_gsva.R

**Purpose:** Run GSVA on gene sets

**Inputs:**
- Filtered expression data from step 01
- `EstroGene_Signatures.xlsx` (Early/Mid/Late gene sets)
- MSigDB pathways (inflammation, immune)

**Logic:**
1. Load gene set signatures
2. Run GSVA on each dataset
3. Combine pathway scores

**Outputs:**
- `outputs/02_gsva_estrogene.RData`
- `outputs/02_gsva_pathways.RData`

### 03_run_mica.R

**Purpose:** Run MICA concordance analysis

**Inputs:**
- GSVA scores from step 02
- Age group labels

**Logic:**
1. Prepare multi-study data matrices (genes/pathways x samples)
2. Run `mica.full()` with n.perm=500
3. Extract gMI+ statistics and p-values
4. Run post-hoc pairwise analysis

**Parameters:**
- `n.perm = 500`
- `p.threshold = 0.05`
- `n.parallel = 10`

**Outputs:**
- `outputs/03_mica_results.RData`
- `outputs/03_mica_summary.csv`

### 04_fig4_estrogen.R

**Purpose:** Generate Figure 4 panels

**Figure 4C: HSD17B Gene Heatmap**
- Genes: ESR1, GREB1, PGR, SAA1, RAB19, KRT37, TRPM8, CYP19A1, HSD17B1/7/12/2/10/14
- Median expression per age group per study
- ComplexHeatmap with row split by study

**Figure 4H: EstroGene Pathway Boxplots**
- Pathways: EstroGene_Early, EstroGene_Mid, EstroGene_Late
- Boxplots faceted by study
- Color by age group

**Outputs:**
- `figures/fig4c_hsd17b_heatmap.pdf`
- `figures/fig4c_hsd17b_heatmap.png`
- `figures/fig4h_estrogen_boxplots.pdf`
- `figures/fig4h_estrogen_boxplots.png`

### 05_fig5_immune.R

**Purpose:** Generate Figure 5 panels

**Figure 5C: Inflammatory Pathway Heatmap**
- Pathways from Neil's edited list (inflammation subset)
- Median GSVA scores per age group per study

**Figure 5E: Immune Cell Boxplots**
- Cell types: CD8 T cells, Cytotoxic lymphocytes, Monocytic lineage, Myeloid DCs
- MCP-counter deconvolution scores
- Boxplots faceted by study

**Outputs:**
- `figures/fig5c_inflammatory_heatmap.pdf`
- `figures/fig5c_inflammatory_heatmap.png`
- `figures/fig5e_immune_boxplots.pdf`
- `figures/fig5e_immune_boxplots.png`

## Code Standards

Per project CLAUDE.md:
- `set.seed(12345)` for reproducibility
- Save figures as PDF + PNG (300 DPI)
- Use `tryCatch()` for error handling
- Export CSV summaries for key results
- Save RData checkpoints

## SLURM Configuration

```bash
#SBATCH --partition=htc
#SBATCH --time=4:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
```

## Required Data Files

See `Jian_MICA/required_data_files.md` for complete list of input files needed from Jian.

## Implementation Plan

1. Create directory structure
2. Implement 01_prep_data.R (once data files obtained)
3. Implement 02_run_gsva.R
4. Implement 03_run_mica.R
5. Implement 04_fig4_estrogen.R
6. Implement 05_fig5_immune.R
7. Create run_all.sbatch
8. Test end-to-end pipeline
9. Validate figures against manuscript
