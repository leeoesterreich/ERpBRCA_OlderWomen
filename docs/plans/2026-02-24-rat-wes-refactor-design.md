# RatWES Refactor to Unified Repo Design

**Date:** 2026-02-24
**Status:** Approved
**Author:** Claude (with user approval)

## Overview

Refactor the Rat Whole Exome Sequencing (WES) analysis from the original Jupyter notebook (`NeilRatWES/Neil_RatWES_Complete.ipynb`) into the unified `ERpBRCA_OlderWomen` repository structure, generating publication-ready figures for the manuscript.

## Requirements

- Create `cleanup` branch from current `main`
- Modify existing scripts in `analysis/03_rat_wes/`
- Generate two manuscript outputs:
  1. **Oncoplot figure** - Mutation landscape heatmap (genes x samples)
  2. **COSMIC signatures plot** - SBS signature decomposition

## Git Workflow

1. Create `cleanup` branch from current `main` in `ERpBRCA_OlderWomen`
2. All refactoring work happens on `cleanup` branch
3. Commit incrementally as each script is updated

## Script Modifications

### `01_parse_vep.py` - No changes needed

Already parses VEP output and filters for HIGH/MODERATE impact variants. Outputs `all_samples_filtered.csv` used by downstream scripts.

### `02_cosmic_signatures.py` - Add figure generation

**Current functionality:**
- Runs SigProfiler COSMIC signature assignment
- Outputs `cosmic_signatures.csv`

**Additions:**
- Generate stacked bar chart (SVG) showing SBS signature counts per sample
- Age group labels on x-axis (Young: 102, 107, 116 / Old: 157, 158, 167)
- Output: `figures/cosmic_signatures.svg`

### `03_generate_oncoplot.py` - Major expansion

**Current functionality:**
- Basic pivot table of Gene x Sample consequences

**Additions:**
- Homolog mapping (rat Ensembl → human Ensembl via pybiomart)
- Gene symbol resolution (human Ensembl → gene symbol)
- Cancer gene filtering (using `brca_genelist.csv`)
- Publication-quality oncoplot figure with:
  - Color-coded consequences (missense, frameshift, etc.)
  - Gene names on y-axis, samples on x-axis
  - Age group labels (Y/O suffix)
  - Legend for consequence types
- Output: `figures/oncoplot.svg`

### `run_analysis.sbatch` - Minor update

- Add `mkdir -p figures`
- Verify conda environment has required packages

## Data Dependencies

### Input files (already exist)

- VEP output: `/bgfs/alee/LO_LAB/General/Lab_Data/20240628_WES_Rat_Neil/results/variant_calling/mutect2/`
- VCF input for SigProfiler: `/bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/RatWES_Mutect2_VCF_Input`

### Reference files to copy into repo

- `brca_genelist.csv` → `analysis/03_rat_wes/data/`

### Sample-to-age mapping

```python
YOUNG_SAMPLES = ['102', '107', '116']  # Young rats
OLD_SAMPLES = ['157', '158', '167']    # Old rats
```

### Expected outputs

```
outputs/
├── all_samples_filtered.csv      # From 01
├── cosmic_signatures.csv         # From 02
├── oncoplot_data.csv            # From 03
├── oncoplot_summary.csv         # From 03
└── homolog_cache.csv            # From 03 (cached mapping)

figures/
├── cosmic_signatures.svg         # From 02
└── oncoplot.svg                  # From 03
```

## Environment & Error Handling

### Conda environment

- Use existing `machine_learning_env`
- Required packages: pandas, matplotlib, seaborn, pybiomart, SigProfilerAssignment

### Error handling

- pybiomart queries: retry logic with chunked queries (200 genes per request)
- SigProfiler: skip re-run if output already exists
- Cache homolog mapping results to avoid repeated API calls

### Figure styling

- Font size: 12-18pt for labels
- Colormap: `tab20` for consequences
- Figure size: oncoplot height scales with gene count
- SVG format for publication quality

## Success Criteria

1. `cleanup` branch created with all changes
2. Pipeline runs successfully via `sbatch run_analysis.sbatch`
3. Both figures generated in `figures/` directory
4. Figures match manuscript style (age group labels, color coding)
