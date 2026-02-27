# Figure 7 Reproducibility Fix Design

**Date:** 2026-02-26
**Status:** Approved
**Author:** Claude (with user approval)

## Overview

Fix Figure 7 panels B/C/D to match the manuscript format exactly, validating that the original analysis code and biostatistics are correct through visual reproducibility.

## Problem

The manuscript validation pipeline flagged Figure 7 panels as FAIL:

| Panel | Manuscript Shows | Current Output | Issue |
|-------|------------------|----------------|-------|
| B/C | Dual side-by-side heatmaps (HALLMARK + BIOCARTA) with raw z-scores | Single difference heatmap (Elderly - Young) | Wrong visualization approach |
| D | CellPhoneDB dot plot (L-R pairs, dot size = p-value, color = expression) | Bar chart fallback | CellChat not installed |

## Goal

Regenerate figures that visually match the manuscript to confirm the underlying code and statistics are sound. This is a reproducibility validation exercise.

## Approach

**Approach 1 (Selected):** Modify existing R scripts + install CellChat

- Minimal code changes
- Reuses existing GSVA calculations
- CellChat produces exact manuscript format via `netVisual_bubble()`

## Design

### Section 1: CellChat Installation

**Environment:** `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging`

**Steps:**
1. Activate conda environment
2. Install CellChat dependencies via conda (if missing): `igraph`, `NMF`, `ggalluvial`
3. Install CellChat from GitHub: `devtools::install_github("jinworks/CellChat")`

**Validation:** `library(CellChat)` loads without error

**Fallback:** If CellChat fails, switch to Python CellPhoneDB (Approach 3)

### Section 2: Panel B/C Heatmap Refactoring

**File:** `analysis/04_human_scrnaseq/15_multicelltype_pathway.R`

**Manuscript Format:**
- Rows: CellType_AgeGroup (e.g., Bcells_Younger, Bcells_Older, Macrophage_Younger, Macrophage_Older)
- Columns: Pathway names
- Two panels: HALLMARK (left), BIOCARTA (right)
- Values: Raw GSVA z-scores (NOT differences)
- Row grouping: Clustered blocks by category (Lymphocyte, Myeloid, Epithelial, Stromal) with visual separation

**Changes:**
1. Keep raw GSVA scores - remove Elderly-Young difference calculation
2. Transpose matrix: cell_type_age (rows) x pathways (columns)
3. Split into two gene set lists: `hallmark_subset`, `biocarta_subset`
4. Row ordering by category with `gaps_row` for visual separation
5. Combine two heatmaps side-by-side using `grid.arrange()` or `cowplot::plot_grid()`

**Output:** `figures/fig7bc_pathway_heatmaps.png`

### Section 3: Panel D CellChat Dot Plot

**File:** `analysis/04_human_scrnaseq/16_cellchat_analysis.R`

**Manuscript Format:**
- Rows: Ligand-receptor pairs (CCL2-CCR2, CXCL9-CXCR3, IL15-IL15RA, etc.)
- Columns: Cell type pairs (sender -> receiver)
- Dot size: -log10(p-value)
- Dot color: Log2 mean expression (yellow-blue scale)

**Changes:**
1. Remove fallback bar chart code
2. Ensure CellChat runs for both age groups
3. Fix `netVisual_bubble()` parameters:
   - `sources.use`: macrophage cell types
   - `targets.use`: T cell/NK cell types
   - `comparison = c(1, 2)` for Young vs Elderly
4. Match dot plot aesthetics to manuscript

**Output:** `figures/fig7d_cellchat_dotplot.png`

### Section 4: Validation & Testing

**Validation approach:**
1. Re-run `04_visual_comparison.py` with updated figure paths
2. Target: Fig. 7_B/C and Fig. 7_D change from FAIL to PASS/SIMILAR

**Statistical verification:**
- Confirm GSVA z-scores use same normalization (log2 CPM)
- Confirm CellChat uses same database (Secreted Signaling)
- Verify cell type annotations match manuscript methods

**Outputs to compare:**
- `multicelltype_pathway_scores.csv` - raw GSVA scores
- `cellchat_*.rds` - communication probabilities

**Success criteria:**
- Visual validation returns PASS or SIMILAR for all Figure 7 panels
- GPT-5 vision confirms "same scientific conclusions"
- Numerical differences within tolerance (rtol=1e-3)

## Files to Modify

| File | Action |
|------|--------|
| `analysis/04_human_scrnaseq/15_multicelltype_pathway.R` | Refactor for dual heatmaps |
| `analysis/04_human_scrnaseq/16_cellchat_analysis.R` | Fix dot plot generation |
| `erp_brca_aging` conda env | Install CellChat |

## Scripts to Run

1. CellChat installation (one-time)
2. `15_multicelltype_pathway.R` - regenerate Panel B/C
3. `16_cellchat_analysis.R` - regenerate Panel D
4. `scripts/validation/04_visual_comparison.py` - validate results
