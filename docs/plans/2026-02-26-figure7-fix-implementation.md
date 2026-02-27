# Figure 7 Reproducibility Fix Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix Figure 7 panels B/C/D to match manuscript format, validating original analysis correctness.

**Architecture:** Modify two R scripts (15_multicelltype_pathway.R, 16_cellchat_analysis.R) to produce manuscript-matching visualizations. Install CellChat for Panel D dot plots.

**Tech Stack:** R, Seurat, GSVA, msigdbr, pheatmap, CellChat, conda

---

## Task 1: Install CellChat in Conda Environment

**Files:**
- Modify: conda environment `erp_brca_aging`

**Step 1: Check if CellChat already installed**

```bash
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging
Rscript -e "library(CellChat); cat('CellChat version:', as.character(packageVersion('CellChat')), '\n')"
```

Expected: Either version number (skip to Task 2) or error "there is no package called 'CellChat'"

**Step 2: Install CellChat dependencies**

```bash
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging
Rscript -e "
if (!require('NMF', quietly=TRUE)) install.packages('NMF', repos='https://cloud.r-project.org')
if (!require('ggalluvial', quietly=TRUE)) install.packages('ggalluvial', repos='https://cloud.r-project.org')
if (!require('ComplexHeatmap', quietly=TRUE)) BiocManager::install('ComplexHeatmap')
"
```

Expected: Packages install or already present

**Step 3: Install CellChat from GitHub**

```bash
Rscript -e "
if (!require('devtools', quietly=TRUE)) install.packages('devtools', repos='https://cloud.r-project.org')
devtools::install_github('jinworks/CellChat')
"
```

Expected: CellChat installs successfully

**Step 4: Verify CellChat installation**

```bash
Rscript -e "library(CellChat); cat('CellChat version:', as.character(packageVersion('CellChat')), '\n')"
```

Expected: Prints version number (e.g., "CellChat version: 2.1.0")

**Step 5: Commit environment notes**

No git commit needed - conda environment is external to repo.

---

## Task 2: Refactor Panel B/C Heatmap Script

**Files:**
- Modify: `analysis/04_human_scrnaseq/15_multicelltype_pathway.R`
- Output: `analysis/04_human_scrnaseq/figures/fig7bc_pathway_heatmaps.png`

**Step 1: Add cowplot/gridExtra to library imports**

In `15_multicelltype_pathway.R`, add after line 23:

```r
  library(gridExtra)
```

**Step 2: Remove difference calculation (lines 150-177)**

Delete or comment out the entire "Step 5: Calculate Elderly - Young difference" section.

**Step 3: Replace Step 6 with dual heatmap generation**

Replace lines 179-226 with:

```r
# -----------------------------------------------------------------------------
# Step 5: Generate Figure 7 B/C - Dual Heatmaps (HALLMARK + BIOCARTA)
# -----------------------------------------------------------------------------
cat("\nStep 5: Generating Figure 7 B/C dual heatmaps...\n")

# Transpose: rows = cell_type_age, columns = pathways
gsva_t <- t(gsva_result)

# Create row annotations for cell type categories
get_category <- function(ct_age) {
  ct <- gsub("_(Elderly|Young)$", "", ct_age)
  case_when(
    grepl("Tcells|NK|Bcells|Plasma", ct) ~ "Lymphocyte",
    grepl("Macro|Mono|DC|Myeloid", ct) ~ "Myeloid",
    grepl("Epithelial|Cancer|Luminal|Basal", ct) ~ "Epithelial",
    grepl("CAF|PVL|Endo|Fibro", ct) ~ "Stromal",
    TRUE ~ "Other"
  )
}

row_categories <- data.frame(
  Category = sapply(rownames(gsva_t), get_category),
  row.names = rownames(gsva_t)
)

# Define row order: group by category, then by cell type (Younger before Older)
category_order <- c("Lymphocyte", "Myeloid", "Epithelial", "Stromal", "Other")
row_order <- rownames(gsva_t)[order(
  match(row_categories$Category, category_order),
  gsub("_(Elderly|Young)$", "", rownames(gsva_t)),
  grepl("_Elderly$", rownames(gsva_t))
)]
gsva_t <- gsva_t[row_order, ]
row_categories <- row_categories[row_order, , drop = FALSE]

# Calculate gaps for category separation
category_counts <- table(row_categories$Category)[category_order]
category_counts <- category_counts[!is.na(category_counts) & category_counts > 0]
gaps_row <- cumsum(category_counts)[-length(category_counts)]

# Split into HALLMARK and BIOCARTA
hallmark_cols <- grep("^HALLMARK_", colnames(gsva_t), value = TRUE)
biocarta_cols <- grep("^BIOCARTA_", colnames(gsva_t), value = TRUE)

gsva_hallmark <- gsva_t[, hallmark_cols, drop = FALSE]
gsva_biocarta <- gsva_t[, biocarta_cols, drop = FALSE]

# Clean column names for display
colnames(gsva_hallmark) <- gsub("^HALLMARK_", "", colnames(gsva_hallmark))
colnames(gsva_hallmark) <- gsub("_", " ", colnames(gsva_hallmark))
colnames(gsva_biocarta) <- gsub("^BIOCARTA_", "", colnames(gsva_biocarta))
colnames(gsva_biocarta) <- gsub("_", " ", colnames(gsva_biocarta))

# Annotation colors
ann_colors <- list(
  Category = c(
    Lymphocyte = "#4DAF4A",
    Myeloid = "#E41A1C",
    Epithelial = "#377EB8",
    Stromal = "#984EA3",
    Other = "grey70"
  )
)

# Color scale (blue-white-red, symmetric around 0)
max_val <- max(abs(gsva_t), na.rm = TRUE)
color_breaks <- seq(-max_val, max_val, length.out = 101)
color_palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Create HALLMARK heatmap
p_hallmark <- pheatmap(
  gsva_hallmark,
  annotation_row = row_categories,
  annotation_colors = ann_colors,
  main = "HALLMARK",
  color = color_palette,
  breaks = color_breaks,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  gaps_row = gaps_row,
  fontsize_row = 8,
  fontsize_col = 9,
  show_rownames = TRUE,
  angle_col = 45,
  silent = TRUE
)

# Create BIOCARTA heatmap
p_biocarta <- pheatmap(
  gsva_biocarta,
  annotation_row = row_categories,
  annotation_colors = ann_colors,
  main = "BIOCARTA",
  color = color_palette,
  breaks = color_breaks,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  gaps_row = gaps_row,
  fontsize_row = 8,
  fontsize_col = 9,
  show_rownames = TRUE,
  angle_col = 45,
  silent = TRUE
)

# Combine side-by-side
png(file.path(figures_dir, "fig7bc_pathway_heatmaps.png"),
    width = 16*300, height = 12*300, res = 300)
grid.arrange(p_hallmark$gtable, p_biocarta$gtable, ncol = 2)
dev.off()
cat("  Saved fig7bc_pathway_heatmaps.png\n")

cat("\n=== Multi-cell-type pathway analysis complete ===\n")
```

**Step 4: Run script to verify it works**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq
sbatch run_figure7_reproduction.sbatch
```

Expected: Job completes, `figures/fig7bc_pathway_heatmaps.png` generated

**Step 5: Commit changes**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
git add analysis/04_human_scrnaseq/15_multicelltype_pathway.R
git commit -m "fix(fig7): refactor heatmap to match manuscript format

- Transpose matrix: rows = cell_type_age, columns = pathways
- Split into HALLMARK and BIOCARTA side-by-side panels
- Use raw GSVA z-scores (not differences)
- Group rows by cell category with visual gaps"
```

---

## Task 3: Fix Panel D CellChat Dot Plot Script

**Files:**
- Modify: `analysis/04_human_scrnaseq/16_cellchat_analysis.R`
- Output: `analysis/04_human_scrnaseq/figures/fig7d_cellchat_dotplot.png`

**Step 1: Remove the fallback bar chart code**

Delete lines 62-106 (the entire `if (!cellchat_available)` block including the `quit()` call).

**Step 2: Fix the netVisual_bubble call (lines 231-241)**

Replace the existing bubble plot section with:

```r
# -----------------------------------------------------------------------------
# Step 5: Generate Figure 7D - Dot Plot
# -----------------------------------------------------------------------------
cat("\nStep 5: Generating Figure 7D...\n")

# Create merged CellChat object for comparison
object.list <- list(Young = cellchat_young, Elderly = cellchat_elderly)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# Identify macrophage and T/NK cell types present in data
all_celltypes <- union(levels(cellchat_young@idents), levels(cellchat_elderly@idents))
macro_types <- grep("Macro|Mono", all_celltypes, value = TRUE)
tcell_types <- grep("Tcells|NK", all_celltypes, value = TRUE)

cat("  Source cell types (Macrophage/Monocyte):", paste(macro_types, collapse = ", "), "\n")
cat("  Target cell types (T/NK cells):", paste(tcell_types, collapse = ", "), "\n")

if (length(macro_types) > 0 && length(tcell_types) > 0) {
  # Generate dot plot
  png(file.path(figures_dir, "fig7d_cellchat_dotplot.png"),
      width = 14*300, height = 12*300, res = 300)

  netVisual_bubble(cellchat_merged,
                   sources.use = macro_types,
                   targets.use = tcell_types,
                   comparison = c(1, 2),
                   angle.x = 45,
                   remove.isolate = TRUE,
                   title.name = "Macrophage/Monocyte -> T/NK Cell Communication")

  dev.off()
  cat("  Saved fig7d_cellchat_dotplot.png\n")
} else {
  cat("  WARNING: Could not find matching cell types for dot plot\n")
}

cat("\n=== CellChat analysis complete ===\n")
```

**Step 3: Run script to verify it works**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq
sbatch run_cellchat_only.sbatch
```

Expected: Job completes, `figures/fig7d_cellchat_dotplot.png` generated with dot plot format

**Step 4: Commit changes**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
git add analysis/04_human_scrnaseq/16_cellchat_analysis.R
git commit -m "fix(fig7): generate proper CellChat dot plot

- Remove fallback bar chart code
- Use netVisual_bubble() for L-R pair visualization
- Show Macrophage -> T/NK cell communication
- Compare Young vs Elderly side-by-side"
```

---

## Task 4: Update Visual Comparison Config

**Files:**
- Modify: `scripts/validation/config.py` (if needed)
- Run: `scripts/validation/04_visual_comparison.py`

**Step 1: Check current figure mapping**

```bash
grep -A5 "Fig.*7" /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/scripts/validation/config.py
```

Expected: See current mapping for Figure 7 panels

**Step 2: Update figure mapping if needed**

If paths need updating, modify config.py to point to new output files:
- `fig7bc_pathway_heatmaps.png` for Panel B/C
- `fig7d_cellchat_dotplot.png` for Panel D

**Step 3: Run visual comparison**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
python scripts/validation/04_visual_comparison.py
```

Expected: Figure 7 panels change from FAIL to PASS or SIMILAR

**Step 4: Review validation report**

```bash
cat validation/report/validation_report.md | grep -A10 "Fig.*7"
```

Expected: All Figure 7 panels show PASS or SIMILAR with notes about matching scientific conclusions

**Step 5: Commit validation updates**

```bash
git add scripts/validation/config.py validation/
git commit -m "fix(validation): update Figure 7 paths and re-run validation

Figure 7 panels B/C/D now match manuscript format.
Visual validation confirms same scientific conclusions."
```

---

## Task 5: Final Verification

**Step 1: View regenerated figures side-by-side with manuscript**

Use the Read tool to visually inspect:
- `analysis/04_human_scrnaseq/figures/fig7bc_pathway_heatmaps.png`
- `analysis/04_human_scrnaseq/figures/fig7d_cellchat_dotplot.png`

Compare against manuscript originals:
- `validation/extracted/figures/main/slide_07_Figure 7/img_02.png`
- `validation/extracted/figures/main/slide_07_Figure 7/img_01.png`

**Step 2: Confirm validation report shows success**

Check that all Figure 7 panels are PASS or SIMILAR in the final validation report.

**Step 3: Final commit if any adjustments needed**

```bash
git add -A
git commit -m "fix(fig7): final adjustments for manuscript match"
```

---

## Summary

| Task | Description | Output |
|------|-------------|--------|
| 1 | Install CellChat | CellChat available in R |
| 2 | Refactor Panel B/C | `fig7bc_pathway_heatmaps.png` |
| 3 | Fix Panel D | `fig7d_cellchat_dotplot.png` |
| 4 | Update validation | Updated `validation_report.md` |
| 5 | Final verification | All panels PASS/SIMILAR |

**Fallback:** If CellChat installation fails, switch to Python CellPhoneDB approach (see design doc for details).
