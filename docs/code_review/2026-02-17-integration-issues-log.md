# Code Integration Issues Log

**Date:** 2026-02-17
**Reviewer:** Claude (AI-assisted biostatistical review)

## Summary

| Category | Issues Found | Fixed |
|----------|-------------|-------|
| P-value correction | 3 | 3 |
| Code bugs | 6 | 6 |
| Missing code | 3 | 3 |
| HPC policy | 2 | 2 |

## Issues Fixed

### Human scRNA-seq (GitHub)

| Issue | File | Fix |
|-------|------|-----|
| No FDR correction on cell fraction tests | 04_cell_fractions.R | Added `p.adjust(method = "BH")` |
| Undefined `SeuratObj_ERpos_Young/Elderly` | 09_cellphonedb_prep.R | Properly subset using `AgeGroup` |
| Typo `CountData_YoungElderlyroWhole` | 09_cellphonedb_prep.R | Fixed variable name |
| Typo `_Macrophage` vs `_MonoMacro` | 03_cell_type_annotation.R | Check `CellTypeMinor` exists before use |
| Missing GSVA code | - | Created 06_run_gsva.R (pseudo-bulk) |
| Missing PROGENy code | - | Created 07_run_progeny.R (pseudo-bulk) |
| Missing WCSEA code | - | Created 08_run_wcsea.R |

### Rat Bulk RNA-seq (Zip)

| Issue | File | Fix |
|-------|------|-----|
| DESeq2 results not FDR-filtered | 03_deseq2.R | Filter by `padj < 0.05` not raw p-value |
| `input_dir1` vs `input_dir` variable | 01_alignment.sh | Changed to `INPUT_DIR` consistently |
| Missing SBATCH mail directives | All .sh files | Added `--mail-type=FAIL --mail-user=` |
| `annot.nkis` inappropriate for rat | 04_pam50_subtyping.R | Create minimal annotation data frame |
| HTSeq missing strand flag | 02_htseq_count.sh | Added `-s $STRAND` parameter |

## Reproducibility Gaps Filled

1. **scRNA GSVA**: Created pseudo-bulk approach aggregating per patient
2. **scRNA PROGENy**: Same pseudo-bulk approach
3. **WCSEA**: Implemented using indepthPathway package, using DE genes as input

## Biostatistical Best Practices Applied

- All multiple comparison tests now use BH-FDR correction
- Significant results filtered at padj < 0.05 (not raw p-value)
- Non-parametric tests used where sample sizes are small (Wilcoxon)
- Gene coverage checks added for cross-species PAM50 subtyping
