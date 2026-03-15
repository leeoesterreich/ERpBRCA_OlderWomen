# ERpBRCA_OlderWomen Project Guidelines

## Code Quality Standards

**All code in this project must be robust, rigorous, and reproducible.**

### Robustness Requirements

1. **Error Handling**: Wrap potentially failing operations in `tryCatch()` (R) or `try/except` (Python)
2. **Graceful Degradation**: If a non-critical step fails, log a warning and continue
3. **Input Validation**: Check that input files exist before processing
4. **Version Compatibility**: Handle API deprecations and package version conflicts explicitly

### Known Compatibility Issues

**Seurat 5.x + patchwork S4 bug:**
- The `+` operator fails when combining `VlnPlot()` with `plot_annotation()` due to S4 method dispatch
- **Fix**: Use `wrap_plots()` to combine plots before adding annotations
- Reference: https://github.com/satijalab/seurat/issues/7653
- See `05_gene_expression_violin.R` for the `safe_vlnplot()` wrapper pattern

**decoupleR 2.x API changes:**
- `dc.get_dorothea()` → `dc.op.dorothea()`
- `dc.decouple()` → `dc.mt.decouple()` or `dc.mt.ulm()`
- Results: `adata.obsm['score_ulm']` not `adata.obsm['ulm_estimate']`

### Reproducibility Requirements

1. **Set seeds**: All scripts must set `set.seed(12345)` (R) or `np.random.seed(12345)` (Python)
   - **Exception**: Section `07_organoid_single_cell` uses `seed=42` / `random_state=42` to preserve consistency with original analysis results. Do not change to 12345.
2. **Fixed ordering**: Use explicit ordering for pathways, genes, cell types - never rely on hash order
3. **Version pinning**: Document package versions in environment files
4. **Deterministic outputs**: Same inputs must produce identical outputs

### Output Standards

1. **Dual format**: Save figures as both PDF (vector) and PNG (300 DPI for validation)
2. **Intermediate files**: Save RDS/pickle checkpoints for long-running analyses
3. **CSV exports**: Save key results as CSV for easy inspection and cross-platform access

## Project Structure

```
analysis/
├── 01_human_bulk_rnaseq/   # Human ER+ bulk RNA-seq (MICA, correlations)
├── 02_rat_snrnaseq/        # Rat mammary tumor snRNA-seq
├── 03_rat_wes/             # Rat whole exome sequencing
├── 04_human_scrnaseq/      # Human ER+ scRNA-seq (Wu et al.)
├── 05_rat_bulk_rnaseq/     # Rat bulk RNA-seq
├── 06_spatial_biopsies/   # Spatial biopsy immune analysis (HCC22-088)
├── 07_organoid_single_cell/ # Organoid scRNA-seq (PDO-296 HSD17B7i)
└── .pipeline_markers/      # Job completion markers

validation/
├── extracted/              # Figures/claims from manuscript
├── comparisons/            # Validation results (numerical, visual, claims)
└── report/                 # Validation reports
```

## Validation Pipeline

Before submission, run the 3-tier validation:
1. **Numerical**: Programmatic comparison of p-values, correlations, sample sizes
2. **Visual**: Azure GPT-5 vision comparison of manuscript vs regenerated figures
3. **Claims**: Text-to-data mapping verification

## Conda Environment

**Primary environment**: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging`

Contains: R 4.3.3, Seurat 5.3.0, GSVA, msigdbr, CellChat, decoupleR, DESeq2

## Duplicate Sample Warning

**`135938ERPRnoNAC` is a duplicate of `1410354A13KECKHE`** - always exclude from analyses.
