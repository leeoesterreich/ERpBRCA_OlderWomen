# Design: Publication-Ready Code & Methods Fixes

**Date:** 2026-03-19
**Goal:** Fix all code issues and methods discrepancies identified in the methods-vs-code audit, making ERpBRCA_OlderWomen a rigorous, reproducible, publication-ready monorepo.

---

## 1. Quick Code Fixes

One-liner changes across 8 files. No behavioral impact on results.

| Issue | File | Change |
|-------|------|--------|
| W3 | `04_human_scrnaseq/15_multicelltype_pathway.R:8` | Header comment "Wilcoxon" → "Welch t-test" |
| W13 | `03_rat_wes/03_generate_oncoplot.py:217` | Fix comment: Young=(157,158,167), Old=(102,107,116) |
| W26 | `01_human_bulk_rnaseq/05_visualize.R` | Remove 3 duplicate `ggsave()` calls (lines 96, 124, 143) |
| W27 | `02_rat_snrnaseq/03b_sctype_annotate.R:85` | Fix `GenesIn` → `GeneIn` variable name |
| W29 | `03_rat_wes/03_generate_oncoplot.py:187` | `.applymap()` → `.map()` |
| W35 | `05_rat_bulk_rnaseq/04_pam50_subtyping.R` | Remove redundant `"NUF2"="Nuf2"` alias |
| W36 | `05_rat_bulk_rnaseq/04_pam50_subtyping.R` | Add PNG 300 DPI output for both heatmaps |
| W20 | `06_spatial_biopsies/02_immune_pathways.py:99` | Remove deprecated `organism='Human'` from `gp.enrichr()` |

Update `code_review_issues.md`: strikethrough + RESOLVED for W3, W13, W20, W26, W27, W29, W35, W36.

---

## 2. Missing Random Seeds

Add seeds to files that lack them per project convention. Section 07 uses `seed=42`; all others use `12345`.

| File | Seed | What to add |
|------|------|-------------|
| `03_rat_wes/01_parse_vep.py` | 12345 | `np.random.seed(12345)` at top |
| `03_rat_wes/02_cosmic_signatures.py` | 12345 | `np.random.seed(12345)` at top |
| `04_human_scrnaseq/16_run_cellphonedb.py` | 42 | Add `debug_seed=42` to `cpdb_statistical_analysis_method.call()` |
| `02_rat_snrnaseq/01_qc_filter.R` | 12345 | Add `seed.use=12345` to per-sample `SCTransform()` in DoubletFinder loop |
| `02_rat_snrnaseq/02_normalize_integrate.R` | 12345 | Add explicit `seed` param to `RunHarmony()` if available, else document reliance on global seed |
| `07_organoid_single_cell/02_qc.py` | 42 | `np.random.seed(42)` at top |
| `07_organoid_single_cell/03_preprocess.py` | 42 | `np.random.seed(42)` + `random_state=42` on PCA/UMAP/Leiden calls |
| `07_organoid_single_cell/06_cell_cycle.py` | 42 | `np.random.seed(42)` at top |
| `07_organoid_single_cell/07_single_cell_pathways.py` | 42 | `np.random.seed(42)` at top |
| `07_organoid_single_cell/08_inhibitor_mechanism_exploration.py` | 42 | `np.random.seed(42)` at top |

Update `code_review_issues.md`: strikethrough + RESOLVED for W14, W15, W21, W28.

---

## 3. Substantive Code Fixes

### 3a. H1 -- Unify curated pathway sets (script 15)

**File:** `04_human_scrnaseq/15_multicelltype_pathway.R`

The FDR correction set (lines ~237-243) and the visualization set (lines ~393-401) must be unified. Use the visualization set (which includes IL2_STAT5, IFN_ALPHA, EMT; excludes COMPLEMENT) for both FDR correction and display.

**Change:** Replace the FDR curated pathway list with the visualization curated list. This means FDR correction will be applied across the biologically curated 25 pathways consistently.

### 3b. W24 -- Duplicate gene removal by highest expression

**File:** `01_human_bulk_rnaseq/01_preprocess.R`

Replace `!duplicated(GeneSymb)` (arbitrary first-occurrence) with selection of the row with highest mean expression across samples. Apply to both TPM and count modes.

Pattern:
```r
# For each duplicated gene symbol, keep the row with highest mean expression
dup_genes <- unique(gene_symbols[duplicated(gene_symbols)])
keep <- rep(TRUE, nrow(expr_matrix))
for (g in dup_genes) {
  idx <- which(gene_symbols == g)
  means <- rowMeans(expr_matrix[idx, , drop = FALSE])
  keep[idx] <- FALSE
  keep[idx[which.max(means)]] <- TRUE
}
```

Log the number of duplicates resolved.

### 3c. M18 -- Add annotation to count-mode heatmap

**File:** `01_human_bulk_rnaseq/01_preprocess.R`

Add `annotation_col` parameter to the count-mode distance heatmap `pheatmap()` call (around line 413-420), matching the TPM-mode correlation heatmap's age range and tissue type annotations.

Update `code_review_issues.md`: strikethrough + RESOLVED for W24. Add new entry for H1 as resolved.

---

## 4. VEP Cache Update + WES Re-run

### 4a. Download cache v104

Download `rattus_norvegicus_vep_104_Rnor_6.0.tar.gz` from Ensembl FTP and install alongside the existing v95 cache.

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/
curl -O https://ftp.ensembl.org/pub/release-104/variation/vep/rattus_norvegicus_vep_104_Rnor_6.0.tar.gz
tar xzf rattus_norvegicus_vep_104_Rnor_6.0.tar.gz
```

### 4b. Update VEP sbatch

**File:** `03_rat_wes/00_run_vep.sbatch`

Change `--cache_version 95` → `--cache_version 104` and update the cache directory path.

### 4c. Re-run WES pipeline

Submit `00_run_vep.sbatch` then `run_analysis.sbatch`. The pipeline is: VEP annotation → parse/filter → COSMIC signatures → oncoplot.

### 4d. Spot-check results

Compare new vs old:
- `cosmic_signatures.csv` -- same dominant signatures?
- `oncoplot_data.csv` -- same genes, same consequence types?
- Visual comparison of oncoplot and signature bar chart

If results are concordant: proceed. If materially different: flag for discussion before continuing.

Update `code_review_issues.md`: strikethrough + RESOLVED for W6.

---

## 5. MICA Integration

### 5a. Adapt Jian_MICA scripts into monorepo

Copy and adapt the working scripts from `Jian_MICA/scripts/` into `analysis/01_human_bulk_rnaseq/mica/`, replacing the stubs. The existing `01_prep_data.R` has real `create_age_groups()` logic (42 lines) that must be reconciled with the Jian version, not blindly overwritten.

Root-level scripts in `Jian_MICA/` (`03_240408_figure_update.R`, `04_EstroGene_MICA[1].R`, `04_EstroGene_MICA.R`) are older versions superseded by `scripts/` -- ignore them.

| Jian_MICA source | Target in mica/ | Adaptation needed |
|-----------------|-----------------|-------------------|
| `scripts/00_setup_env.sh` | `00_setup_env.sh` | Adapt env path to monorepo convention |
| `scripts/00_unzip_data.sh` | `00_unzip_data.sh` | Update data paths |
| `scripts/01_estrogene_gsva.R` | `01_prep_data.R` (reconcile with existing) | Merge Jian's implementation with existing age group defs; update data paths |
| `scripts/02_mica_analysis.R` | `03_run_mica.R` (replace stub) | Update input/output paths |
| `scripts/03_figure_generation.R` | `04_fig4_estrogen.R` + `05_fig5_immune.R` (replace stubs) | Split into figure-specific scripts; merge helper functions from existing stubs if needed |
| `scripts/04_extract_claims.py` | `06_extract_claims.py` | Adapt paths |
| `scripts/05_statistical_audit.R` | `07_statistical_audit.R` | Adapt paths |
| `scripts/06_compile_report.py` | `08_compile_report.py` | Adapt paths |

### 5b. Copy required data

Copy data files from `Jian_MICA/data/` into `analysis/01_human_bulk_rnaseq/mica/data/input/`. Total size: ~1.3 GB (dominated by `ilc_clean_data.RData` at 1.2 GB).

Directory structure to copy:
```
data/
├── clinical/       # TCGA/METABRIC clinical info
├── database/       # MICA database files
├── gmt/            # Gene set GMT files
├── ilc_clean_data.RData  # Pre-compiled expression data (1.2 GB)
└── raw_data/       # TCGA/METABRIC/SCAN-B raw matrices
```

Reference `Jian_MICA/required_data_files.md` for expected inputs.

### 5c. Update sbatch

Adapt `Jian_MICA/sbatch_*.sh` (4 scripts: gsva, mica, figures, audit) into a unified `mica/run_all.sbatch` that uses the `erp_brca_aging` environment and follows monorepo conventions (mail directives, completion markers).

### 5d. Verify execution

Run the integrated MICA pipeline and confirm outputs match `Jian_MICA/results/` and `Jian_MICA/figures/`.

### 5e. Add `set.seed(12345)` to any MICA scripts missing it

Update `code_review_issues.md`: strikethrough + RESOLVED for C3, W10, W25.

---

## 6. Methods Documentation Corrections

Update all 7 section-level methods files + main text + code_review_issues.md.

### 6a. Section 01 (`01_human_bulk_rnaseq_methods.md`)

- Remove "Log2(TPM)" from PROGENy input description (H6); state it uses VST matrix
- Add numeric age group boundaries for primary analysis (M20/W22)
- Document MICA integration (no longer "planned" -- now executed)
- Note duplicate gene handling change (now highest expression, not first occurrence)

### 6b. Section 02 (`02_rat_snrnaseq_methods.md`)

- Fix DE cell filter: "fewer than 10 cells per sample, requiring >=2 samples per group" (H2)
- State explicitly that DE/DA uses marker-scoring annotation at resolution 0.4, not scType (H3/W11)
- Fix conda env name: `erp_snrnaseq` not `erp_brca_aging` (M1)
- Add gene filter before DESeq2: >=10 counts in >=2 samples (M13)
- Document propeller re-correction rationale (M12)
- Fix DoubletFinder PCA: document that 50 PCs are computed (default), with dims 1:30 used for paramSweep and doubletFinder (M4)

### 6c. Section 03 (`03_rat_wes_methods.md`)

- Update VEP cache version from 95 to 104 (after re-run)
- Fix conda env name: `aging_wes` not `erp_brca_aging` (M2)
- Note that original analysis used VEP v95; results regenerated with v114.2 + cache v104

### 6d. Section 04 (`04_human_scrnaseq_methods.md`)

- Verify R version and correct table (H5) -- check `erp_brca_aging` env
- Document that script 15 uses min_patients=2 vs script 14 uses min_patients=3 (M10)
- Fix "top50" → "top40" or note actual count (M11)
- Clarify script 13 produces DEG bar chart, not L-R interaction analysis (M17)
- Document unified curated pathway set after H1 fix
- Note script 09 GSEA uses Wilcoxon score ranking and 3 libraries (M14-15)

### 6e. Section 05 (`05_rat_bulk_rnaseq_methods.md`)

- Fix validation matching description: hardcoded lookup for HTSeq, gsub for PAM50 (M19)

### 6f. Section 06 (`06_spatial_biopsies_methods.md`)

- Change "combined" heatmap → "separate enrichment and prerank heatmaps" (M16)
- Add CITEgeist version/parameters or cite the tool paper (V2)
- Add Carleton et al. full citation with year, journal, DOI (V3)

### 6g. Section 07 (`07_organoid_scrnaseq_methods.md`)

- Add Cell Ranger version 9.0.1 (M5)
- Document script 09 GSEA uses Wilcoxon score ranking and 3 libraries (M14-15)
- Add reference build details (Ensembl 109 FASTA, GENCODE v44 GTF) (M6)
- Add brief note that pseudobulk DE has n=1 replicate per condition, acknowledging statistical limitation (V6)

### 6h. Main text (`00_main_text_methods.md`)

- Update section 01 summary to reflect MICA integration
- Update section 03 to reflect VEP cache v104

### 6i. Code review issues (`code_review/code_review_issues.md`)

- Remove W18 (wrong project -- ERT Visium duplicate sample, confirmed not in bulk RNA-seq input data)
- Strikethrough + RESOLVED for all issues fixed in batches 1-5, including: W3, W6, W11, W13, W14, W15, W20, W21, W24, W25, W26, W27, W28, W29, W35, W36, C3, W10
- Add new resolved entry for H1 (curated pathway set unification)
- Keep remaining open issues (W4, W7, W8, W12, W16, W17, W19, W22, W23, W30-34) as accepted limitations

---

## Execution Order

1. **Batch 1** (quick code fixes) -- no dependencies, safe to parallelize
2. **Batch 2** (seeds) -- no dependencies, safe to parallelize with batch 1
3. **Batch 3** (substantive fixes) -- no dependencies on 1/2
4. **Batch 4** (VEP cache + re-run) -- independent, can run in parallel with 1-3; requires HPC job
5. **Batch 5** (MICA integration) -- independent of 1-4; requires file copying + HPC job to verify
6. **Batch 6** (methods docs) -- depends on all code fixes being complete; update after batches 1-5

Batches 1-3 can be done as parallel subagent work (different files, no conflicts). Batch 4 requires a SLURM job. Batch 5 requires reading Jian_MICA scripts and adapting. Batch 6 is the final pass.

---

## Out of Scope (Accepted Limitations)

These remain open but are documented, not blocking:

- W4: GSVA dense matrix OOM risk (operational, not a bug)
- W7: `aging_wes` env spec untracked (add yml to repo if desired)
- W8/W16: BioMart not pinned to archive (cache files mitigate)
- W12: Fixed 5% doublet rate (standard practice)
- W17: STAR genome index provenance (institutional resource)
- W19: Package versions unpinned (capture sessionInfo() separately)
- W22/W23: Age group boundary consistency, TPM fallback filename (document)
- W30-34: Hardcoded paths, SLURM dependencies (institutional conventions)

### Low-Priority Audit Items (Not Addressed)

These are cosmetic or minor documentation gaps that do not affect reproducibility or correctness:

- L1-L2: Undocumented CSV export step and density plot in section 01 QC
- L3: TMB, K-means, CNVKit absent from refactored WES pipeline (intentionally dropped)
- L4: SigProfilerMatrixGenerator genome install prerequisite undocumented
- L5-L8: Figure formatting inconsistencies in section 06 (DPI, dual format, font sizes)
- L9-L10: Undocumented fulvestrant calibration and convex hull analyses in section 07
- L11: Normalized counts in DESeq2 results CSV not mentioned in methods
- L12: GSVA heatmap colormap described as blue-white-red but is blue-white-orange-red
- V1: E1UpRegGene "407 genes" may differ after cleaning
- V4-V5: UMAP "default parameters" and reference genome vagueness (addressed partially by M6)
- M7: Methods tables list "(not pinned)" for package versions (addressed with W19 if sessionInfo captured)
