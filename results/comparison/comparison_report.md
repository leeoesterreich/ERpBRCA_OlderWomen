# ERpBRCA_OlderWomen: Biostatistical Corrections Summary

Generated: 2026-02-06

## Overview

This report documents the biostatistical corrections applied during code review
and validates that the core results are reproducible.

## Human Bulk RNA-seq Analysis

### Corrections Applied:
1. **FDR Correction**: Added Benjamini-Hochberg FDR adjustment for multiple testing
2. **Random Seeds**: Added for reproducibility (seed = 12345)

### Results Comparison:
| Metric | Original | Corrected |
|--------|----------|-----------|
| Total tests | 50 | 50 |
| Significant (p < 0.05) | 22 | 22 |
| Significant (FDR < 0.05) | N/A | 20 |

### Reproducibility Status:
**PASS**: All 22 nominally significant correlations are identical between original and corrected analyses.

After FDR correction, 2 tests lost significance:
- HSD17B7 × Estrogen (p=0.026, FDR=0.061)
- HSD17B7 × HALLMARK_ESTROGEN_RESPONSE_EARLY (p=0.036, FDR=0.081)

## Rat snRNA-seq Analysis

### Corrections Applied:
1. **Differential Expression**: Added Young vs Aged comparison per cell type (WAS MISSING)
2. **Differential Abundance**: Added propeller testing for cell type proportions (WAS MISSING)
3. **DoubletFinder**: Added doublet detection
4. **FDR Correction**: Applied BH correction to all tests
5. **Random Seeds**: Added for reproducibility
6. **JoinLayers()**: Added for Seurat v5 compatibility

### New Results Generated:

**Differential Expression (Young vs Aged):**

| Cell Type | Total DE Genes | FDR < 0.05 | Up in Aged | Down in Aged |
|-----------|----------------|------------|------------|--------------|
| Luminal | 3,482 | 3,335 | 978 | 2,357 |
| Basal | 3,802 | 2,682 | 1,811 | 871 |
| Fibroblasts | 4,747 | 1,631 | 1,479 | 152 |
| Endothelial | 4,308 | 1,112 | 645 | 467 |
| **Total** | **16,339** | **8,760** | **4,913** | **3,847** |

**Differential Abundance:**
- Cell types tested: 4
- Significant at FDR < 0.05: 0
- Note: Limited statistical power with n=3 per group

## Technical Fixes Applied

1. **sys.frame(1)$ofile error**: Replaced with `commandArgs()` for Rscript compatibility
2. **future.globals.maxSize**: Increased to 4GB for large cell populations
3. **Cluster name 'g' prefix**: Fixed mapping for Seurat v5 numeric cluster names
4. **tryCatch scope**: Fixed variable scoping in error handlers

## Conclusion

The biostatistical corrections ensure:
- ✓ Proper multiple testing correction (FDR)
- ✓ Complete statistical analyses (DE/DA for rat data)
- ✓ Reproducibility through random seeds
- ✓ Seurat v5 compatibility

**Core correlation values from the human bulk RNA-seq analysis are fully reproduced**,
demonstrating that the underlying data processing is consistent. The FDR correction
appropriately reduces the number of significant findings from 22 to 20.

The rat snRNA-seq analysis now includes comprehensive DE and DA testing that was
previously missing, with 8,760 FDR-significant differentially expressed genes
identified across 4 cell types.
