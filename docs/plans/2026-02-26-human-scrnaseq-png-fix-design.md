# Human scRNA-seq PNG Figure Fix Design

**Date:** 2026-02-26
**Status:** Approved
**Author:** Claude (with user approval)

## Overview

Fix the human scRNA-seq pipeline to generate PNG figures and update the validation pipeline to find and compare them against manuscript figures.

## Problem

The validation report shows only 4 visual comparisons were performed, with just one human scRNA-seq figure (umap_celltypes.png → Fig. 7_B). The R scripts generate PDFs but the validation pipeline requires PNGs. Additional human scRNA-seq figures (fraction_boxplot, gsva_heatmap, progeny_heatmap) aren't being validated.

## Solution

Two-pronged approach:
1. **Primary:** Modify R scripts to save PNG directly to `figures/` directory
2. **Fallback:** Add PDF-to-PNG conversion step in sbatch script for any missed files
3. **Validation:** Update `FIGURE_MAPPING` in `04_visual_comparison.py` to map all human scRNA-seq figures

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    HUMAN SCRNASEQ PNG FIX                       │
├─────────────────────────────────────────────────────────────────┤
│  ┌──────────────────┐    ┌──────────────────┐                  │
│  │ R Scripts        │    │ Conversion       │                  │
│  │ (Primary)        │    │ Script           │                  │
│  │                  │    │ (Fallback)       │                  │
│  │ - 03_cell_type   │    │                  │                  │
│  │ - 04_cell_frac   │    │ PDF → PNG        │                  │
│  │ - 06_gsva        │    │ for any missed   │                  │
│  │ - 07_progeny     │    │                  │                  │
│  └────────┬─────────┘    └────────┬─────────┘                  │
│           │                       │                             │
│           ▼                       ▼                             │
│  ┌─────────────────────────────────────────┐                   │
│  │     figures/                             │                   │
│  │     ├── umap_celltypes.png              │                   │
│  │     ├── fraction_boxplot.png            │                   │
│  │     ├── gsva_heatmap.png                │                   │
│  │     └── progeny_heatmap.png             │                   │
│  └─────────────────────────────────────────┘                   │
│                         │                                       │
│                         ▼                                       │
│  ┌─────────────────────────────────────────┐                   │
│  │ Validation Pipeline                      │                   │
│  │ 04_visual_comparison.py                  │                   │
│  │                                          │                   │
│  │ FIGURE_MAPPING updated for all panels   │                   │
│  └─────────────────────────────────────────┘                   │
└─────────────────────────────────────────────────────────────────┘
```

## Files to Modify

| File | Change |
|------|--------|
| `analysis/04_human_scrnaseq/03_cell_type_annotation.R` | Add ggsave for UMAP |
| `analysis/04_human_scrnaseq/04_cell_fractions.R` | Add ggsave for boxplot |
| `analysis/04_human_scrnaseq/06_run_gsva.R` | Add png() device for heatmap |
| `analysis/04_human_scrnaseq/07_run_progeny.R` | Add png() device for heatmap |
| `analysis/04_human_scrnaseq/run_analysis.sbatch` | Add figures dir + conversion step |
| `scripts/validation/04_visual_comparison.py` | Update FIGURE_MAPPING |

## Implementation Details

### R Script Changes

Each script needs:
1. Create figures directory: `figures_dir <- "figures"; dir.create(figures_dir, showWarnings = FALSE)`
2. Add PNG save after existing PDF block

**03_cell_type_annotation.R** (after line 92):
```r
ggsave(file.path(figures_dir, "umap_celltypes.png"), p1 + p2,
       width = 14, height = 6, dpi = 300, bg = "white")
```

**04_cell_fractions.R**:
```r
ggsave(file.path(figures_dir, "fraction_boxplot.png"), p,
       width = 12, height = 6, dpi = 300, bg = "white")
```

**06_run_gsva.R** (pheatmap needs png() device):
```r
png(file.path(figures_dir, "gsva_heatmap.png"), width = 10*300, height = 6*300, res = 300)
print(p)
dev.off()
```

**07_run_progeny.R**:
```r
png(file.path(figures_dir, "progeny_heatmap.png"), width = 10*300, height = 8*300, res = 300)
print(p)
dev.off()
```

### Fallback Conversion (run_analysis.sbatch)

Add at end before completion marker:
```bash
echo "=== Step 11: Ensure PNG figures exist ==="
for pdf in outputs/*.pdf; do
    [ -f "$pdf" ] || continue
    base=$(basename "$pdf" .pdf)
    png="figures/${base}.png"
    if [ ! -f "$png" ]; then
        echo "Converting $pdf -> $png"
        env -i PATH=/usr/bin:/bin pdftoppm -png -r 300 -singlefile "$pdf" "figures/${base}" 2>/dev/null \
            || env -i PATH=/usr/bin:/bin convert -density 300 "$pdf" "$png" 2>/dev/null \
            || echo "WARNING: Could not convert $pdf"
    fi
done
```

Note: `env -i` avoids conda library conflicts with system tools.

### Validation Mapping Update

Update `FIGURE_MAPPING` in `scripts/validation/04_visual_comparison.py`:

```python
"slide_07_Figure 7": {
    "figure_id": "Fig. 7",
    "analysis": "04_human_scrnaseq",
    "source": "main",
    "panels": {
        "B_celltypes": {
            "regenerated": "figures/umap_celltypes.png",
            "manuscript_hint": ["img_01", "img_02"],
        },
        "C_pathways": {
            "regenerated": "figures/gsva_heatmap.png",
            "manuscript_hint": ["img_06", "img_07"],
        },
        "D_fractions": {
            "regenerated": "figures/fraction_boxplot.png",
            "manuscript_hint": ["img_08", "img_09"],
        },
        "E_progeny": {
            "regenerated": "figures/progeny_heatmap.png",
            "manuscript_hint": ["img_10", "img_11"],
        },
    }
},
```

## Testing

1. Run the modified sbatch: `sbatch run_analysis.sbatch`
2. Verify PNG files exist in `figures/` directory
3. Run validation: `python scripts/validation/04_visual_comparison.py`
4. Check that all human scRNA-seq figures appear in `visual_comparison.json`

## Notes

- PNG output at 300 DPI for publication quality
- `bg = "white"` in ggsave prevents transparent backgrounds
- pheatmap objects require png() device (ggsave doesn't work reliably with them)
- Fallback conversion uses pdftoppm (faster) with ImageMagick convert as backup
