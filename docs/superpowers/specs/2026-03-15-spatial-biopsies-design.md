# Design: Section 06 — Spatial Biopsies

**Date:** 2026-03-15
**Status:** Approved

## Context

Incorporate spatial biopsy immune analysis (scripts 12 and 13) from the CITEgeist/Neil
analysis into the ERpBRCA_OlderWomen project as a new analysis section. These are
HCC22-088 samples from a separate Carleton et al. paper, used here as supporting
spatial evidence for the inflammation/macrophage story. The other paper will be cited
in the manuscript.

## Approach

**Minimal copy** — copy scripts and dependencies into the new section, update paths
and env references, no logic changes.

## Structure

```
analysis/06_spatial_biopsies/
├── 01_immune_secretion.py      # Spatial cytokine plots (from CITEgeist script 12)
├── 02_immune_pathways.py       # Enrichment + GSEA (from CITEgeist script 13)
├── utils.py                    # log_memory_usage only (extracted, not full copy)
├── figure_config.py            # Standardized figure params (copied from CITEgeist)
├── data/
│   └── biopsy_adatas.pkl       # Symlink to CITEgeist data
├── run_analysis.sbatch         # Combined sbatch runner
├── run_analysis.sh             # Shell runner script
├── logs/
└── figures/
    ├── immune_secretion/
    └── immune_pathways/
        ├── enrichment/
        ├── prerank/
        └── dotplots/
```

## Changes from Originals

1. **Renumbered**: 12 → 01, 13 → 02 (per-section convention)
2. **Log filenames**: `logs/12_immune_secretion.log` → `logs/01_immune_secretion.log`,
   `logs/13_immune_pathways.log` → `logs/02_immune_pathways.log`
3. **Conda env**: sbatch points to `erp_brca_aging` instead of `CITEgeist_env_neil`
4. **sbatch**: Single `run_analysis.sbatch` running both scripts sequentially, with
   `cd` into section directory, mail directives per project standards. Supports
   optional `--figures-only` flag passed through to script 02.
5. **Data**: Symlink from `data/biopsy_adatas.pkl` to CITEgeist source
6. **utils.py**: Extract only `log_memory_usage()` (~5 lines, needs `psutil`, `os`,
   `logging`) rather than copying the full 600-line file with unneeded dependencies
7. **figure_config.py docstring**: Update "CITEgeist" reference to current project
8. **GSEA seed**: Change `seed=42` to `seed=12345` in `gp.prerank()` call to match
   CLAUDE.md reproducibility standard
9. **No other analysis logic changes**: Cytokine lists, sample mappings, pathway
   filters all kept identical

## Dependencies

- Python packages needed in `erp_brca_aging`: scanpy, gseapy, seaborn, psutil
  (verify/install)
- Data: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/data/biopsy_adatas.pkl`
- Font: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Arial.ttf` (already on cluster)

## Known Issues (Accepted)

From Codex review — accepted as-is for minimal copy approach:

- **obsm shape mismatch** (script 13 line 431): likely dead code path for this dataset;
  no obsm keys in biopsy_adatas match macrophage keywords
- **Scanpy show_ prefix workaround** (script 12): brittle but functional
- **No input file validation**: violates CLAUDE.md but scripts already work
- **Hardcoded Arial path**: fine on this cluster
- **Hardcoded 6-sample order**: correct for this dataset
- **Figure format**: Script 12 saves PNG only, script 13 saves PNG+SVG (not PDF).
  CLAUDE.md asks for PDF+PNG but these are spatial plots where SVG/PNG is sufficient.
  Accepted deviation.

## CLAUDE.md Update

Add `06_spatial_biopsies` to the project structure section in CLAUDE.md.
