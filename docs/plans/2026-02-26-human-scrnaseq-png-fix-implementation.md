# Human scRNA-seq PNG Figure Fix Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add PNG figure generation to human scRNA-seq pipeline and update validation mapping.

**Architecture:** Modify 4 R scripts to save PNG to `figures/` directory, add fallback PDF-to-PNG conversion to sbatch, update validation `FIGURE_MAPPING` to find all new figures.

**Tech Stack:** R (ggsave, png device), Bash (pdftoppm, ImageMagick convert), Python (validation script)

---

## Task 1: Add figures directory setup to 03_cell_type_annotation.R

**Files:**
- Modify: `analysis/04_human_scrnaseq/03_cell_type_annotation.R:33-34`

**Step 1: Add figures_dir variable after output_dir**

After line 33 (`output_dir <- file.path(script_dir, "outputs")`), add:

```r
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
```

**Step 2: Add ggsave call after dev.off()**

After line 92 (`dev.off()`), add:

```r
# Save PNG for validation pipeline
ggsave(file.path(figures_dir, "umap_celltypes.png"), p1 + p2,
       width = 14, height = 6, dpi = 300, bg = "white")
```

**Step 3: Verify by visual inspection**

Run: `grep -n "figures_dir\|ggsave" analysis/04_human_scrnaseq/03_cell_type_annotation.R`
Expected: Lines showing figures_dir creation and ggsave call

**Step 4: Commit**

```bash
git add analysis/04_human_scrnaseq/03_cell_type_annotation.R
git commit -m "feat(human-scrnaseq): add PNG output to cell type annotation"
```

---

## Task 2: Add PNG output to 04_cell_fractions.R

**Files:**
- Modify: `analysis/04_human_scrnaseq/04_cell_fractions.R:33-34,143`

**Step 1: Add figures_dir variable after output_dir**

After line 34 (`output_dir <- file.path(script_dir, "outputs")`), add:

```r
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
```

**Step 2: Add ggsave call after dev.off()**

After line 143 (`dev.off()`), add:

```r
# Save PNG for validation pipeline
ggsave(file.path(figures_dir, "fraction_boxplot.png"), p,
       width = 12, height = 6, dpi = 300, bg = "white")
```

**Step 3: Verify by visual inspection**

Run: `grep -n "figures_dir\|ggsave" analysis/04_human_scrnaseq/04_cell_fractions.R`
Expected: Lines showing figures_dir creation and ggsave call

**Step 4: Commit**

```bash
git add analysis/04_human_scrnaseq/04_cell_fractions.R
git commit -m "feat(human-scrnaseq): add PNG output to cell fractions"
```

---

## Task 3: Add PNG output to 06_run_gsva.R

**Files:**
- Modify: `analysis/04_human_scrnaseq/06_run_gsva.R:35-36,148`

**Step 1: Add figures_dir variable after output_dir**

After line 35 (`output_dir <- file.path(script_dir, "outputs")`), add:

```r
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
```

**Step 2: Add png() device call after dev.off()**

After line 148 (`dev.off()`), add:

```r
# Save PNG for validation pipeline (pheatmap needs explicit device)
png(file.path(figures_dir, "gsva_heatmap.png"), width = 10*300, height = 6*300, res = 300)
pheatmap(
  gsva_result,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  main = "GSVA Estrogen Pathway Scores (Pseudo-bulk)",
  scale = "row",
  show_colnames = TRUE
)
dev.off()
```

**Step 3: Verify by visual inspection**

Run: `grep -n "figures_dir\|png(" analysis/04_human_scrnaseq/06_run_gsva.R`
Expected: Lines showing figures_dir creation and png() device call

**Step 4: Commit**

```bash
git add analysis/04_human_scrnaseq/06_run_gsva.R
git commit -m "feat(human-scrnaseq): add PNG output to GSVA heatmap"
```

---

## Task 4: Add PNG output to 07_run_progeny.R

**Files:**
- Modify: `analysis/04_human_scrnaseq/07_run_progeny.R:32-33,111`

**Step 1: Add figures_dir variable after output_dir**

After line 32 (`output_dir <- file.path(script_dir, "outputs")`), add:

```r
figures_dir <- file.path(script_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
```

**Step 2: Add png() device call after dev.off()**

After line 111 (`dev.off()`), add:

```r
# Save PNG for validation pipeline (pheatmap needs explicit device)
png(file.path(figures_dir, "progeny_heatmap.png"), width = 10*300, height = 8*300, res = 300)
pheatmap(
  pathway_activity,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  main = "PROGENy Pathway Activity (Pseudo-bulk)",
  scale = "column",
  cluster_cols = FALSE
)
dev.off()
```

**Step 3: Verify by visual inspection**

Run: `grep -n "figures_dir\|png(" analysis/04_human_scrnaseq/07_run_progeny.R`
Expected: Lines showing figures_dir creation and png() device call

**Step 4: Commit**

```bash
git add analysis/04_human_scrnaseq/07_run_progeny.R
git commit -m "feat(human-scrnaseq): add PNG output to PROGENy heatmap"
```

---

## Task 5: Add fallback PDF-to-PNG conversion to run_analysis.sbatch

**Files:**
- Modify: `analysis/04_human_scrnaseq/run_analysis.sbatch:26,115-121`

**Step 1: Add figures to mkdir command**

Change line 26 from:
```bash
mkdir -p logs outputs outputs/preprocessing
```
to:
```bash
mkdir -p logs outputs outputs/preprocessing figures
```

**Step 2: Add conversion step before completion marker**

Before line 117 (`MARKER_DIR=...`), add:

```bash
echo "=== Step 11: Ensure PNG figures exist ==="
# Convert any PDFs that don't have PNG counterparts (fallback)
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

**Step 3: Verify by visual inspection**

Run: `grep -n "figures\|pdftoppm\|convert" analysis/04_human_scrnaseq/run_analysis.sbatch`
Expected: Lines showing figures mkdir, pdftoppm and convert calls

**Step 4: Commit**

```bash
git add analysis/04_human_scrnaseq/run_analysis.sbatch
git commit -m "feat(human-scrnaseq): add PDF-to-PNG fallback conversion"
```

---

## Task 6: Update validation FIGURE_MAPPING

**Files:**
- Modify: `scripts/validation/04_visual_comparison.py:157-168`

**Step 1: Expand slide_07_Figure 7 panels**

Replace lines 157-168 (the entire `"slide_07_Figure 7"` entry) with:

```python
        "slide_07_Figure 7": {
            "figure_id": "Fig. 7",
            "analysis": "04_human_scrnaseq",
            "source": "main",
            "panels": {
                # Panel B: Cell type UMAP showing transcriptomic differences
                "B_celltypes": {
                    "regenerated": "figures/umap_celltypes.png",
                    "manuscript_hint": ["img_01", "img_02"],
                },
                # Panel C: DEG/pathway analysis - GSVA heatmap
                "C_pathways": {
                    "regenerated": "figures/gsva_heatmap.png",
                    "manuscript_hint": ["img_06", "img_07"],
                },
                # Panel D: Cell fractions comparison
                "D_fractions": {
                    "regenerated": "figures/fraction_boxplot.png",
                    "manuscript_hint": ["img_08", "img_09"],
                },
                # PROGENy pathway activity
                "E_progeny": {
                    "regenerated": "figures/progeny_heatmap.png",
                    "manuscript_hint": ["img_10", "img_11"],
                },
            }
        },
```

**Step 2: Verify by visual inspection**

Run: `grep -A 30 "slide_07_Figure 7" scripts/validation/04_visual_comparison.py`
Expected: All 4 panels (B_celltypes, C_pathways, D_fractions, E_progeny) listed

**Step 3: Commit**

```bash
git add scripts/validation/04_visual_comparison.py
git commit -m "feat(validation): add human scRNA-seq figure mappings"
```

---

## Task 7: Integration test - run pipeline and validation

**Step 1: Submit the human scRNA-seq pipeline**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq
sbatch run_analysis.sbatch
```

Note job ID for monitoring.

**Step 2: Wait for completion and verify figures exist**

After job completes:
```bash
ls -la figures/*.png
```

Expected: 4 PNG files (umap_celltypes.png, fraction_boxplot.png, gsva_heatmap.png, progeny_heatmap.png)

**Step 3: Run visual validation**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
python scripts/validation/04_visual_comparison.py --dry-run
```

Expected: Shows all human scRNA-seq figure pairs found

**Step 4: Run full validation**

```bash
python scripts/validation/04_visual_comparison.py
```

Expected: visual_comparison.json updated with human scRNA-seq results

**Step 5: Regenerate validation report**

```bash
python scripts/validation/06_generate_report.py
```

**Step 6: Commit validation results (if successful)**

```bash
git add validation/comparisons/visual/visual_comparison.json validation/report/validation_report.md
git commit -m "docs: update validation report with human scRNA-seq figures"
```

---

## Summary

| Task | File | Change |
|------|------|--------|
| 1 | 03_cell_type_annotation.R | Add figures_dir + ggsave |
| 2 | 04_cell_fractions.R | Add figures_dir + ggsave |
| 3 | 06_run_gsva.R | Add figures_dir + png() device |
| 4 | 07_run_progeny.R | Add figures_dir + png() device |
| 5 | run_analysis.sbatch | Add figures mkdir + PDF conversion |
| 6 | 04_visual_comparison.py | Expand FIGURE_MAPPING |
| 7 | Integration test | Run pipeline + validation |

Total: 7 tasks, ~6 file modifications, ~30 minutes implementation time
