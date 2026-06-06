# Externally-Sourced Code (Dev Repo Only — Not for Public Repo)

These generators originally lived outside the ERpBRCA_OlderWomen tree. They are
kept here for reproducibility audit; the public repo doesn't carry them because
they are upstream code from sibling projects (NeilRatWES, CITEgeistNeilAnalysis).

## Bundle source map

| Panel | Bundle source SVG | Ported here | Original |
|---|---|---|---|
| 2A | `bundle/eps_for_neil/2A.eps` | n/a — public repo has the same script under analysis/03_rat_wes/ | `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/output_plots/oncoplot_BRCA.svg` |
| 7I | sender_top_pathways | `analysis/06_spatial_biopsies/05_commot_analysis.py` + `commot_lib/utils.py:plot_single_group_pathway` | `CITEgeistNeilAnalysis/CITEgeist/analysis/05_commot_analysis.py` + `utils.py` |
| 8I | mean_colocalization_heatmap_combined_macrophages | `analysis/06_spatial_biopsies/04_colocalization.py` + `figure_config_citegeist.py` | `CITEgeistNeilAnalysis/CITEgeist/analysis/09_colocalization_analysis.py` |

## Reverse-engineered helpers (also in public repo, but the lineage stays here)

- `analysis/02_rat_snrnaseq/10_multicelltype_pathway.R` — rewrite for Fig 2E rat HALLMARK GSVA. Sanghoon's original generator was never in either repo.
- `analysis/04_human_scrnaseq/15c_fig7d_elderly_minus_young.R` — Fig 7D Elderly−Young diff heatmap. Original output never committed.
- `analysis/06_spatial_biopsies/03_fig7g_immunosuppressive_macrophages.py` — Fig 7G per-patient panels. Replaces the broken empty-SVG output from `01_immune_secretion.py`.
- `analysis/07_organoid_single_cell/11_crop_fig5de.py` — Fig 5D / 5E single-subpanel re-renders.

## Superseded

- `analysis/04_human_scrnaseq/15b_fig2e_pathway_heatmap.R` — first 2E attempt using human scRNA-seq data; superseded once we noticed the accepted figure legend explicitly says rat snRNA-seq.

## Build driver

- `scripts/build_eps_bundle.sh` — the original session driver. Hardcoded internal OneDrive + lab-share paths; uses `REGENERATE` / `RASTER_FALLBACK` flags that are pre-final iteration state. Kept here for archival.
