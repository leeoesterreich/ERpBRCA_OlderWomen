# Repository Unification and Biostatistical Corrections Design

**Date:** 2026-02-05
**Author:** Code Review Session
**Status:** Approved

## Overview

Unify four separate GitHub repositories into a single, reproducible analysis repository for the paper "Chronic inflammation and hormone disposition promote a tumor-permissive locale for breast cancer in older women" (Carleton et al.). Implement biostatistical corrections identified during code review while preserving the ability to compare original and corrected results.

## Repositories to Unify

| Repository | Content | GEO Accession |
|------------|---------|---------------|
| ERpBRCA_OlderWomen | Main overview (becomes unified home) | - |
| BRCA_OlderWomen_BulkRNAseq | Human bulk RNA-seq analysis | GSE276755 |
| BRCA_Rat_YoungerOlder_snRNAseq | Rat snRNA-seq analysis | GSE276758 |
| RatAgingWES | Rat whole exome sequencing | GSE276759 |

**Excluded:** MICA (external R package by jianzou75, to be incorporated later if needed)

## Key Decisions

- **Code format:** Standalone `.R` and `.sh` scripts (extracted from README markdown)
- **Data handling:** Relative paths with GEO download scripts
- **Figure organization:** Generate in analysis folders, map to manuscript structure
- **Dependencies:** Conda environments for all R and Python analyses
- **Biostatistical approach:** Fix and compare (run both original and corrected analyses)

---

## Repository Structure

```
ERpBRCA_OlderWomen/
├── README.md                    # Overview, abstract, citation info
├── environment.yml              # Master conda environment
├── config.yaml                  # User-configurable paths (if needed)
│
├── data/
│   ├── download_all.sh          # Master download script
│   ├── README.md                # Data descriptions, GEO accessions
│   ├── human_bulk_rnaseq/       # GSE276755
│   │   ├── download.sh
│   │   ├── raw/
│   │   └── external/
│   ├── rat_snrnaseq/            # GSE276758
│   │   ├── download.sh
│   │   ├── raw/
│   │   └── metadata/
│   └── rat_wes/                 # GSE276759
│       ├── download.sh
│       └── raw/
│
├── analysis/
│   ├── 01_human_bulk_rnaseq/
│   │   ├── environment.yml
│   │   ├── 00_download_data.sh
│   │   ├── 01_preprocess.R
│   │   ├── 02_run_gsva.R
│   │   ├── 03_run_progeny.R
│   │   ├── 04_correlations.R
│   │   ├── 05_visualize.R
│   │   ├── run_analysis.sh
│   │   └── outputs/
│   │
│   ├── 02_rat_snrnaseq/
│   │   ├── environment.yml
│   │   ├── 00_download_data.sh
│   │   ├── 01_qc_filter.R
│   │   ├── 02_normalize_integrate.R
│   │   ├── 03_cluster_annotate.R
│   │   ├── 04_subset_myeloid.R
│   │   ├── 05_subset_nkt.R
│   │   ├── 06_cell_fractions.R
│   │   ├── 07_visualize.R
│   │   ├── 08_differential_expression.R
│   │   ├── 09_differential_abundance.R
│   │   ├── 10_statistical_summary.R
│   │   ├── run_analysis.sh
│   │   └── outputs/
│   │
│   └── 03_rat_wes/
│       ├── environment.yml
│       ├── 00_download_data.sh
│       ├── 01_run_sarek.sh
│       ├── 02_run_vep.sh
│       ├── 03_analysis.ipynb
│       ├── run_analysis.sh
│       └── outputs/
│
├── figures/
│   ├── by_analysis/
│   │   ├── human_bulk_rnaseq/
│   │   ├── rat_snrnaseq/
│   │   └── rat_wes/
│   └── manuscript/
│       ├── Fig1/
│       ├── Fig2/
│       ├── Fig3/
│       ├── FigS1/
│       └── FigS2/
│
├── results/
│   ├── original/
│   │   ├── human_bulk_rnaseq/
│   │   └── rat_snrnaseq/
│   ├── corrected/
│   │   ├── human_bulk_rnaseq/
│   │   └── rat_snrnaseq/
│   └── comparison/
│       ├── correlation_comparison.csv
│       ├── significant_findings_summary.md
│       ├── comparison_report.html
│       ├── comparison_report.pdf
│       └── figures/
│
├── scripts/
│   ├── run_all.sh
│   ├── compare_results.R
│   └── organize_figures.R
│
├── envs/
│   ├── bulk_rnaseq.yml
│   ├── snrnaseq.yml
│   └── wes.yml
│
└── docs/
    ├── methods.md
    └── plans/
```

---

## Biostatistical Fixes

### Human Bulk RNA-seq

| Issue | Current | Corrected |
|-------|---------|-----------|
| Multiple testing | None | Benjamini-Hochberg FDR on all 88 correlation tests |
| DESeq2 design | `design = ~1` | `design = ~ AgeRange + Group` |
| Random seeds | None | `set.seed(12345)` before PROGENy |
| Gene filtering | `LOC\\d` regex | `^LOC\\d+` to catch all |
| RAB19 removal | Unexplained | Document rationale or remove from both |
| QC plots | Missing | Add PCA, library size, MA plots |
| Effect sizes | No CIs | Add 95% confidence intervals to correlations |

### Rat snRNA-seq

| Issue | Current | Corrected |
|-------|---------|-----------|
| Differential expression | Not performed | `FindMarkers()` Young vs Aged per cell type, Wilcoxon, BH-adjusted |
| Differential abundance | Not performed | `propeller` or `speckle` for cell type proportion testing |
| Doublet detection | Commented out | Enable DoubletFinder |
| Random seeds | None | Seeds before UMAP, Harmony, FindClusters |
| Hardcoded paths | OneDrive paths | Relative `data/` paths |
| Cell type validation | Manual only | Add SingleR or scType reference validation |
| Sample size | Not discussed | Document n=3/group limitation, report effect sizes |
| Violin plots | No statistics | Add statistical annotations (wilcox.test with BH correction) |

---

## Comparison Framework

### Output Structure

```
results/comparison/
├── correlation_comparison.csv       # Side-by-side p vs FDR
├── significant_findings_summary.md  # What remains significant
├── comparison_report.html           # Visual report
├── comparison_report.pdf            # PDF version for sharing
└── figures/
    ├── bubbleplot_sidebyside.pdf
    ├── significance_change_heatmap.pdf
    └── effect_size_comparison.pdf
```

### Comparison Report Contents

1. Executive summary (1 page) - key findings preserved or changed
2. Methods - what corrections were applied
3. Bulk RNA-seq comparison - correlation tables, bubble plots
4. snRNA-seq comparison - clustering concordance, new DE/DA results
5. Appendix - full statistical tables

---

## Conda Environment

### Master Environment (`environment.yml`)

```yaml
name: erp_brca_aging
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  # R and core packages
  - r-base=4.4.1
  - r-tidyverse
  - r-data.table
  - r-ggplot2
  - r-patchwork
  - r-rmarkdown
  - r-pagedown

  # Bioconductor
  - bioconductor-deseq2
  - bioconductor-gsva
  - bioconductor-progeny
  - bioconductor-seurat
  - bioconductor-harmony
  - bioconductor-singler
  - bioconductor-glmgampoi
  - bioconductor-dittoseq

  # Python
  - python=3.10
  - pandas
  - numpy
  - jupyter
  - matplotlib
  - seaborn

  # Tools
  - pandoc
  - wget
  - curl
```

---

## Execution Flow

### Master Script (`scripts/run_all.sh`)

```bash
#!/bin/bash
set -euo pipefail

echo "=== ERpBRCA OlderWomen Analysis Pipeline ==="

# Step 1: Download data
echo "[1/6] Downloading data from GEO..."
./data/download_all.sh

# Step 2: Run analyses
echo "[2/6] Running human bulk RNA-seq analysis..."
cd analysis/01_human_bulk_rnaseq && ./run_analysis.sh && cd ../..

echo "[3/6] Running rat snRNA-seq analysis..."
cd analysis/02_rat_snrnaseq && ./run_analysis.sh && cd ../..

echo "[4/6] Running rat WES analysis..."
cd analysis/03_rat_wes && ./run_analysis.sh && cd ../..

# Step 3: Compare original vs corrected
echo "[5/6] Generating comparison report..."
Rscript scripts/compare_results.R

# Step 4: Organize figures
echo "[6/6] Organizing figures..."
Rscript scripts/organize_figures.R

echo "=== Pipeline complete ==="
echo "Results: results/comparison/comparison_report.pdf"
echo "Figures: figures/manuscript/"
```

---

## Implementation Plan

### Phase 1: Repository Setup
1. Restructure ERpBRCA_OlderWomen with new directory layout
2. Create conda environment files
3. Write data download scripts
4. Update main README with overview and instructions

### Phase 2: Extract and Refactor Code
5. Extract bulk RNA-seq code from README to standalone `.R` scripts
6. Extract snRNA-seq code from README to standalone `.R` scripts
7. Clean up WES notebook and scripts
8. Add proper headers, random seeds, and logging to all scripts

### Phase 3: Implement Biostatistical Fixes
9. Bulk RNA-seq: Add FDR correction, fix DESeq2 design
10. snRNA-seq: Enable doublet detection, add DE/DA testing
11. Add QC visualizations to both analyses
12. Create comparison framework scripts

### Phase 4: Validation and Figures
13. Run both original and corrected analyses
14. Generate comparison report
15. Organize figures for manuscript
16. Verify all figures can be regenerated cleanly

### Phase 5: Documentation
17. Write detailed methods for supplement
18. Document any result changes for revision letter
19. Final README polish and usage instructions

---

## Success Criteria

- [ ] All analyses run from a single `./scripts/run_all.sh` command
- [ ] Data downloads automatically from GEO
- [ ] Random seeds ensure identical results on re-run
- [ ] FDR-corrected results available alongside original
- [ ] Comparison report clearly shows what changed
- [ ] All manuscript figures generated in `figures/manuscript/`
- [ ] Conda environment installs without errors
- [ ] README provides clear setup and usage instructions
