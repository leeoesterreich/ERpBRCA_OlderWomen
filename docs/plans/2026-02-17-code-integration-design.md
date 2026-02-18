# Code Integration Design: Human scRNA-seq & Rat Bulk RNA-seq

**Date:** 2026-02-17
**Status:** Approved
**Author:** Claude (with user approval)

## Overview

Integrate two new code sources into the consolidated ERpBRCA_OlderWomen repository, including biostatistical review and creation of missing analysis code for reproducibility.

## Sources

### 1. Human scRNA-seq (GitHub)
- **Source:** https://github.com/leeoesterreich/BRCA_Human_YoungerOlder_scRNAseq
- **Content:** Seurat-based analysis of Wu et al. 2021 ER+ breast cancer scRNA-seq data
- **Samples:** 10 patients (3 Young, 4 MidAge, 3 Elderly), ~31,000 cells

### 2. Rat Bulk RNA-seq (Zip)
- **Source:** `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/Re_ Bulk RNA-seq for rat tumors.zip`
- **Content:** STAR alignment, HTSeq counting, DESeq2 DE analysis, PAM50 subtyping
- **Scripts:** Alignment.sh, HTSeq.sh, DEseq2.R, Genefu.R

## Target Structure

```
analysis/
├── 04_human_scrnaseq/
│   ├── 01_load_subset_data.R
│   ├── 02_harmony_integrate.R
│   ├── 03_cell_type_annotation.R
│   ├── 04_cell_fractions.R
│   ├── 05_gene_expression_violin.R
│   ├── 06_run_gsva.R           # NEW - adapted from bulk
│   ├── 07_run_progeny.R        # NEW - adapted from bulk
│   ├── 08_run_wcsea.R          # NEW - from indepthPathway
│   ├── 09_cellphonedb_prep.R
│   ├── run_analysis.sbatch
│   └── run_analysis.sh
│
└── 05_rat_bulk_rnaseq/
    ├── 01_alignment.sh
    ├── 02_htseq_count.sh
    ├── 03_deseq2.R
    ├── 04_pam50_subtyping.R
    ├── run_analysis.sbatch
    └── run_analysis.sh
```

## Review Checklist (Light Review)

### Biostatistical Issues
- [ ] Correct statistical tests for data type
- [ ] Multiple testing correction applied (BH-FDR)
- [ ] Sample size appropriate for test
- [ ] Correct contrast specification in DE analysis

### Code Bugs
- [ ] Variable name typos/mismatches
- [ ] Hardcoded paths parameterized
- [ ] Missing library loads
- [ ] Indexing errors

### P-value Specific
- [ ] FDR correction on multi-comparison tests
- [ ] Correct p.adjust method
- [ ] No uncorrected p-value cherry-picking

## Issues Identified

### Human scRNA-seq (GitHub)
| Issue | Location | Type |
|-------|----------|------|
| No FDR correction on cell fraction comparisons | Step 6 boxplots | P-value |
| Spearman correlations not corrected | Violin loop | P-value |
| Typo: `_Macrophage` vs `_MonoMacro` | Step 3 | Bug |
| Undefined: `SeuratObj_ERpos_Young`/`_Elderly` | Step 6 | Bug |
| Typo: `CountData_YoungElderlyroWhole` | Step 6 | Bug |

### Rat Bulk RNA-seq (Zip)
| Issue | Location | Type |
|-------|----------|------|
| DESeq2 results not FDR-filtered | DEseq2.R:26 | P-value |
| `input_dir1` vs `input_dir` mismatch | Alignment.sh:49,54 | Bug |
| Missing `#SBATCH --mail-*` directives | Both .sh files | HPC policy |
| `annot.nkis` may not match rat data | Genefu.R:73 | Statistical |
| HTSeq missing strand flag | HTSeq.sh:41 | Potential bug |

### Reproducibility Gaps
| Issue | Source | Severity |
|-------|--------|----------|
| scRNA GSVA code missing | Human scRNA-seq | High |
| scRNA PROGENy code missing | Human scRNA-seq | High |
| WCSEA code missing | Human scRNA-seq | Medium |

## Approach

**Sequential Review-Then-Integrate:**

1. **Phase 1: Review** - Check existing code for bugs and p-value issues
2. **Phase 2: Create** - Write missing GSVA/PROGENy/WCSEA code for scRNA
3. **Phase 3: Convert** - Transform to numbered scripts, apply fixes
4. **Phase 4: Verify** - Check paths, dependencies, corrections applied

## Deliverables

- `analysis/04_human_scrnaseq/` - 9+ numbered scripts
- `analysis/05_rat_bulk_rnaseq/` - 4 numbered scripts (2 shell, 2 R)
- Issues log documenting fixes applied
- Updated README.md with new analysis directories
- Updated environment.yml if new dependencies needed

## Dependencies

### Existing (in environment.yml)
- Seurat, harmony, ggplot2, DESeq2, GSVA, progeny

### May Need to Add
- indepthPathway (for WCSEA)
- genefu (for PAM50)

## Notes

- Shell scripts for alignment/counting stay as .sh (SLURM jobs)
- All R scripts get consistent header with library loads
- Match existing repo patterns for run_analysis.sh/sbatch
