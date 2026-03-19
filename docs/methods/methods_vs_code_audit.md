# Methods vs. Code Audit Report

**Date:** 2026-03-19
**Purpose:** Pre-publication verification that all methods documentation accurately reflects the actual code.
**Scope:** All 7 section-level methods files in `docs/methods/` cross-referenced against `analysis/` code.

---

## Executive Summary

All 7 methods sections were reviewed line-by-line against their corresponding analysis code. The methods are **generally accurate and detailed**, with most parameter values, function calls, and thresholds matching precisely. However, **23 discrepancies** were identified that need correction before publication, including 6 high-priority issues that could affect reproducibility claims or reviewer confidence.

### Counts by Section

| Section | Accurate | Discrepancies | Missing from Methods | Missing from Code | Vague |
|---------|----------|---------------|---------------------|-------------------|-------|
| 01 Human Bulk RNA-seq | 25+ | 6 | 9 | 3 | 7 |
| 02 Rat snRNA-seq | 30+ | 9 | 9 | 2 | 6 |
| 03 Rat WES | 20+ | ~5 | 7 | 0 | 5 |
| 04 Human scRNA-seq | 35+ | ~10 | 9 | 1 | 6 |
| 05 Rat Bulk RNA-seq | 25+ | 4 | 7 | 0 | 3 |
| 06 Spatial Biopsies | 20+ | 5 | 8 | 1 | 4 |
| 07 Organoid scRNA-seq | 50+ | 4 | 12 | 0 | 5 |

---

## HIGH-PRIORITY ISSUES (Fix Before Submission)

These discrepancies could be flagged by reviewers or affect reproducibility claims.

### H1. Section 04 -- Two different curated pathway sets for FDR vs. visualization
**File:** `15_multicelltype_pathway.R`
**Problem:** The FDR correction set (lines 237-243) differs from the visualization set (lines 393-401). The FDR set includes `HALLMARK_COMPLEMENT` but omits `IL2_STAT5_SIGNALING`, `INTERFERON_ALPHA_RESPONSE`, and `EPITHELIAL_MESENCHYMAL_TRANSITION`. The methods only documents the visualization set, but statistical conclusions depend on the FDR set.
**Fix:** Document both sets explicitly, or unify them.

### H2. Section 02 -- DE minimum cell threshold described incorrectly
**File:** `08_differential_expression.R`
**Problem:** Methods says "Cell types with fewer than 10 cells in either group were excluded." Code actually filters per-sample (10 cells per sample, then requires >=2 samples per group). This is a fundamentally different statistical design.
**Fix:** Correct the methods to describe the per-sample filter and minimum samples per group.

### H3. Section 02 -- Unclear which annotation feeds DE/DA
**Files:** `03_cluster_annotate.R`, `03b_sctype_annotate.R`, `03c_sctype_cluster.R`, `08_differential_expression.R`
**Problem:** Three annotation strategies are described, but the methods never states which one provides cell type labels for DE/DA. The main pipeline (resolution 0.4) feeds `seurat_annotated.rds` into DE/DA, while scType (resolution 1.5) produces separate objects. Additionally, `CellTypeMacroTcell_RatsnRNAseq` is identical to `CellTypeByMarker_RatsnRNAseq` in the primary pipeline, making the "two levels of granularity" claim in the DA section misleading.
**Fix:** Explicitly state that DE/DA uses the marker-scoring annotation at resolution 0.4.

### H4. Section 03 -- VEP version upgrade undocumented
**Files:** `RatAgingWES/sbatch_vep.sh` (v95) vs `analysis/03_rat_wes/00_run_vep.sbatch` (v114.2)
**Problem:** The original analysis used VEP v95; the refactored pipeline upgraded to v114.2. If published results came from v95, the methods should say v95. If results were regenerated with v114.2, annotation differences should be acknowledged.
**Fix:** Clarify which VEP version produced the published results and note any annotation differences.

### H5. Section 04 -- R version mismatch
**Problem:** Methods table says R 4.4.1; CLAUDE.md says the conda environment has R 4.3.3. This needs to be verified and corrected.
**Fix:** Check `R --version` in the `erp_brca_aging` environment and update the methods accordingly.

### H6. Section 01 -- PROGENy input description self-contradictory
**File:** `01_human_bulk_rnaseq_methods.md`
**Problem:** Line 83-84 says PROGENy uses "Log2(TPM) expression matrix"; line 92 says it uses "the same VST-normalized matrix." Code loads `vst_normalized_matrix.rds`. The two statements contradict each other.
**Fix:** Remove the incorrect TPM reference. State that PROGENy uses the VST matrix (which in TPM-only mode contains log2-TPM).

---

## MEDIUM-PRIORITY ISSUES (Should Fix)

### Environment / Reproducibility

| # | Section | Issue |
|---|---------|-------|
| M1 | 02 | Conda env name wrong: methods says `erp_brca_aging`, sbatch activates `erp_snrnaseq` |
| M2 | 03 | Conda env name wrong: methods says `erp_brca_aging`, sbatch activates `aging_wes` |
| M3 | 02 | DoubletFinder per-sample SCTransform missing `seed.use` (reproducibility claim violation) |
| M4 | 02 | DoubletFinder per-sample PCA computes 50 PCs (default), not 30 as methods implies |
| M5 | 07 | Cell Ranger version (9.0.1) never mentioned in methods |
| M6 | 07 | Reference build details (Ensembl 109 FASTA, GENCODE v44 GTF, biotype filtering, PAR_Y exclusion) not documented |
| M7 | All | Most package versions listed as "(not pinned)" -- need actual version numbers for publication |
| M8 | 07 | Scripts 02, 03, 06, 07, 08 have no random seed despite project convention requiring seed=42 |

### Statistical Methods

| # | Section | Issue |
|---|---------|-------|
| M9 | 04 | Script 15 header says "Wilcoxon" but code runs Welch t-test (W3 in existing issues -- still open) |
| M10 | 04 | Script 15 uses `min_patients_per_group = 2` but script 14 uses `MIN_PATIENTS_PER_GROUP = 3`. Methods doesn't distinguish |
| M11 | 04 | "top50" mode in script 17 actually selects 40 pairs, not 50 |
| M12 | 02 | Propeller BH re-correction: code applies `p.adjust()` on raw p-values, overriding propeller's built-in FDR. Methods should explain why |
| M13 | 02 | Gene-count filter before DESeq2 (>=10 counts in >=2 samples) is undocumented |
| M14 | 07 | Script 09 GSEA uses Wilcoxon `score` for ranking, not `log2FC * -log10(p)` as described for script 05 |
| M15 | 07 | Script 09 GSEA runs against 3 libraries (Hallmark + KEGG + Reactome), not just Hallmarks as implied |

### Factual Errors in Methods Text

| # | Section | Issue |
|---|---------|-------|
| M16 | 06 | Methods claims enrichment + prerank results "combined" into one heatmap; code produces two separate heatmaps |
| M17 | 04 | Script 13 does not perform CellPhoneDB/CellChat analysis as methods implies -- it only makes a bar chart of communication DEGs |
| M18 | 01 | Count-mode distance heatmap is NOT annotated by age/tissue despite methods claiming it is |
| M19 | 05 | Validation matching: methods says "fuzzy name matching" but HTSeq check uses hardcoded lookup table; only PAM50 uses gsub |
| M20 | 01 | Age group boundaries for primary analysis never numerically defined (only MICA sub-analysis has explicit cutoffs) |

---

## LOW-PRIORITY ISSUES (Nice to Fix)

### Missing from Methods (undocumented code steps)

| # | Section | Issue |
|---|---------|-------|
| L1 | 01 | CSV export step in sbatch (correlation_results.csv, mica_results.csv, hsd17b7_correlations.csv) |
| L2 | 01 | Expression distribution density plot in QC |
| L3 | 03 | TMB analysis, K-means clustering, and CNVKit from original notebook absent from refactored pipeline |
| L4 | 03 | SigProfilerMatrixGenerator genome install (`genInstall.install('rn6')`) prerequisite undocumented |
| L5 | 06 | `01_immune_secretion.py` does not use `figure_config.py` -- figures use scanpy defaults, not the 14pt/300 DPI standard |
| L6 | 06 | Font sizes in summary heatmap (9pt y-axis, 11pt x-axis) are below the 14pt minimum claimed in methods |
| L7 | 06 | Spatial plot DPI is likely 150 (scanpy default), not 300 as claimed |
| L8 | 06 | Dual PNG+SVG format not applied to spatial cytokine plots or enrichment barplots |
| L9 | 07 | Fulvestrant calibration and dose-response analyses from script 08 not described |
| L10 | 07 | Convex hull area and sub-clustering from script 10 not described |
| L11 | 05 | Normalized counts merged into DESeq2 results CSV not mentioned |
| L12 | 04 | GSVA heatmap uses blue-white-orange-red, not blue-white-red as methods says |

### Vague Claims Needing Specifics

| # | Section | Issue |
|---|---------|-------|
| V1 | 01 | E1UpRegGene says "407 genes" but actual count after hyphen truncation + deduplication may differ |
| V2 | 06 | No CITEgeist version/parameters/reference provided |
| V3 | 06 | No Carleton et al. citation (year, journal, DOI) |
| V4 | 07 | Reference genome build described as "derived from GRCh38 2024-A" -- insufficient for reproduction |
| V5 | 07 | UMAP "default parameters" relies on knowing scanpy version defaults |
| V6 | 07 | Pseudobulk DE with n=1 replicate per condition has statistical limitations not acknowledged |

---

## Cross-File Consistency Check: Main Text vs. Supplementary

The `00_main_text_methods.md` summaries were compared against the detailed supplementary sections:

| Main Text Claim | Status |
|-----------------|--------|
| Section 01: "top 100 footprint genes per pathway with 1,000 permutations" | Matches code and detailed methods |
| Section 02: "DoubletFinder (5% expected rate, per-sample pK optimization)" | Matches |
| Section 02: "resolution 0.4" | Matches primary pipeline; scType uses 1.5 (not in main text, which is fine) |
| Section 04: "Welch t-tests on per-patient pseudo-bulk scores" | Matches code (but script header says Wilcoxon) |
| Section 04: "CellChat (Secreted Signaling database, truncated mean method)" | Matches (`type = "triMean"`) |
| Section 05: "fastp v0.23.4, Q>=20, length>=36 bp" | Matches |
| Section 06: "top 200 expressed genes" | Matches |
| Section 07: "filter-probes=false" | Matches |
| Section 07: "PyDESeq2, FDR < 0.05, |log2FC| > 0.5" | Matches |

No discrepancies found between main text and supplementary methods.

---

## Relationship to Existing Code Review Issues

Cross-referencing with `docs/code_review/code_review_issues.md`:

| Existing Issue | Status in This Audit |
|----------------|---------------------|
| C1 (DEG column schema) | RESOLVED - methods updated |
| C2 (Multi-cell-type loads macrophage data) | RESOLVED - methods updated |
| C3 (MICA stubs) | Confirmed - methods correctly marks as "planned" |
| W3 (Script 15 Wilcoxon vs Welch) | STILL OPEN - methods says Welch (correct), code header says Wilcoxon (wrong) |
| W6 (VEP cache v95 with software v114.2) | STILL OPEN - see H4 above |
| W7 (SigProfiler not in environment.yml) | Confirmed |
| W8 (BioMart not version-pinned) | Confirmed |
| W11 (Annotation ambiguity) | STILL OPEN - see H3 above |
| W14 (Missing seeds in Python scripts) | Confirmed for sections 03, 04 |
| W21 (Missing seeds in section 07) | Confirmed |

---

## Recommended Action Items

**Before submission:**
1. Fix H1-H6 (high-priority discrepancies)
2. Fix M1-M2 (wrong environment names)
3. Fix M9 (script 15 header comment)
4. Fix M16-M20 (factual errors in methods text)
5. Verify R version (H5) and update methods table
6. Add Cell Ranger version to section 07 methods (M5)

**Nice to have:**
7. Pin all package versions with `sessionInfo()` / `pip freeze` captures
8. Add missing random seeds to section 07 scripts (M8)
9. Document gene counts after cleaning (V1)
10. Add CITEgeist and Carleton et al. citation details (V2-V3)
