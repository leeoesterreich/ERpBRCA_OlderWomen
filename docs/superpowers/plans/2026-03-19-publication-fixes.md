# Publication-Ready Code & Methods Fixes — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix all code issues and methods discrepancies from the methods-vs-code audit, making ERpBRCA_OlderWomen a rigorous, reproducible, publication-ready monorepo.

**Architecture:** Six independent batches of fixes: quick code fixes, random seeds, substantive code changes, VEP cache update + WES re-run, MICA integration, and methods documentation corrections. Batches 1-3 touch different files and can run in parallel. Batch 4 requires an HPC job. Batch 5 requires file copying and adaptation. Batch 6 (methods docs) depends on all prior batches.

**Tech Stack:** R 4.3.3, Python 3.10, Seurat 5.x, scanpy, DESeq2, VEP v114.2, SLURM

**Spec:** `docs/superpowers/specs/2026-03-19-publication-fixes-design.md`
**Audit:** `docs/methods/methods_vs_code_audit.md`
**Issues tracker:** `docs/code_review/code_review_issues.md`

**Base path:** `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen`

---

## Task 1: Quick Code Fixes (Batch 1)

Eight one-liner changes across 6 files. No behavioral impact.

**Files:**
- Modify: `analysis/04_human_scrnaseq/15_multicelltype_pathway.R:8`
- Modify: `analysis/03_rat_wes/03_generate_oncoplot.py:217`
- Modify: `analysis/01_human_bulk_rnaseq/05_visualize.R:96,124,143`
- Modify: `analysis/02_rat_snrnaseq/03b_sctype_annotate.R:85`
- Modify: `analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R`
- Modify: `analysis/06_spatial_biopsies/02_immune_pathways.py:99`

- [ ] **Step 1: Fix W3 — script 15 header comment**

In `analysis/04_human_scrnaseq/15_multicelltype_pathway.R`, line 8: change "Wilcoxon" to "Welch t-test" in the header comment.

- [ ] **Step 2: Fix W13 — oncoplot comment**

In `analysis/03_rat_wes/03_generate_oncoplot.py`, line 217: fix the comment to read `Young (157, 158, 167) = '_Y', Old (102, 107, 116) = '_O'`.

- [ ] **Step 3: Fix W26 — remove duplicate ggsave calls**

In `analysis/01_human_bulk_rnaseq/05_visualize.R`, remove the duplicate `ggsave()` calls at lines 96, 124, and 143 (keep lines 94, 122, 141 respectively).

- [ ] **Step 4: Fix W27 — scType variable typo**

In `analysis/02_rat_snrnaseq/03b_sctype_annotate.R`, line 85: change `GenesIn` to `GeneIn`.

- [ ] **Step 5: Fix W29 — deprecated pandas API**

In `analysis/03_rat_wes/03_generate_oncoplot.py`, line 187: change `.applymap(` to `.map(`.

- [ ] **Step 6: Fix W35 — redundant PAM50 alias**

In `analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R`, find the line `"NUF2" = "Nuf2"` in the ortholog mapping and remove it (keep `"CDCA1" = "Nuf2"`).

- [ ] **Step 7: Fix W36 — add PNG output for PAM50 heatmaps**

In `analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R`, after each `pheatmap()` call that saves to PDF, add PNG output:

For `pam50_heatmap.pdf` (12x15 inches):
```r
png(file.path(figures_dir, "pam50_heatmap.png"), width = 12*300, height = 15*300, res = 300)
# ... same pheatmap call ...
dev.off()
```

For `pam50_probabilities_heatmap.pdf` (12x8 inches):
```r
png(file.path(figures_dir, "pam50_probabilities_heatmap.png"), width = 12*300, height = 8*300, res = 300)
# ... same pheatmap call ...
dev.off()
```

- [ ] **Step 8: Fix W20 — remove deprecated gseapy parameter**

In `analysis/06_spatial_biopsies/02_immune_pathways.py`, line 99: remove `organism='Human'` from the `gp.enrichr()` call.

- [ ] **Step 9: Commit**

```bash
git add analysis/04_human_scrnaseq/15_multicelltype_pathway.R \
        analysis/03_rat_wes/03_generate_oncoplot.py \
        analysis/01_human_bulk_rnaseq/05_visualize.R \
        analysis/02_rat_snrnaseq/03b_sctype_annotate.R \
        analysis/05_rat_bulk_rnaseq/04_pam50_subtyping.R \
        analysis/06_spatial_biopsies/02_immune_pathways.py
git commit -m "fix: quick code fixes for publication review (W3,W13,W20,W26,W27,W29,W35,W36)"
```

---

## Task 2: Missing Random Seeds (Batch 2)

Add reproducibility seeds to 10 files. Section 07 uses `seed=42`; all others use `12345`.

**Files:**
- Modify: `analysis/03_rat_wes/01_parse_vep.py`
- Modify: `analysis/03_rat_wes/02_cosmic_signatures.py`
- Modify: `analysis/04_human_scrnaseq/16_run_cellphonedb.py:91-100`
- Modify: `analysis/02_rat_snrnaseq/01_qc_filter.R:121-126`
- Modify: `analysis/02_rat_snrnaseq/02_normalize_integrate.R:137-143`
- Modify: `analysis/07_organoid_single_cell/02_qc.py`
- Modify: `analysis/07_organoid_single_cell/03_preprocess.py:236,338,359`
- Modify: `analysis/07_organoid_single_cell/06_cell_cycle.py`
- Modify: `analysis/07_organoid_single_cell/07_single_cell_pathways.py`
- Modify: `analysis/07_organoid_single_cell/08_inhibitor_mechanism_exploration.py`

- [ ] **Step 1: Add seeds to WES Python scripts**

In `analysis/03_rat_wes/01_parse_vep.py`, add after imports:
```python
import numpy as np
np.random.seed(12345)
```

In `analysis/03_rat_wes/02_cosmic_signatures.py`, add after imports (numpy is already imported):
```python
np.random.seed(12345)
```

- [ ] **Step 2: Add debug_seed to CellPhoneDB v5**

In `analysis/04_human_scrnaseq/16_run_cellphonedb.py`, lines 91-100, add `debug_seed=42` to the `cpdb_statistical_analysis_method.call()`:
```python
        cpdb_results = cpdb_statistical_analysis_method.call(
            cpdb_file_path=str(CPDB_FILE),
            meta_file_path=str(filtered_meta_path),
            counts_file_path=str(counts_file),
            counts_data='gene_name',
            output_path=str(output_dir),
            threshold=0.1,
            iterations=1000,
            threads=4,
            debug_seed=42,
        )
```

- [ ] **Step 3: Add seed to DoubletFinder SCTransform**

In `analysis/02_rat_snrnaseq/01_qc_filter.R`, lines 121-126, add `seed.use = 12345`:
```r
  seurat_obj <- SCTransform(
    seurat_obj,
    method = "glmGamPoi",
    vars.to.regress = "percent.mt",
    seed.use = 12345,
    verbose = FALSE
  )
```

- [ ] **Step 4: Add seed to RunHarmony**

In `analysis/02_rat_snrnaseq/02_normalize_integrate.R`, lines 137-143. Check if `RunHarmony()` accepts a `seed` parameter. If yes, add it:
```r
seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  assay.use = "SCT",
  plot_convergence = FALSE,
  seed = 12345
)
```
If `RunHarmony()` does not accept `seed`, add a comment: `# Harmony relies on global set.seed(12345) at script top`.

- [ ] **Step 5: Add seeds to section 07 Python scripts**

For each of the following files, add `import numpy as np; np.random.seed(42)` (or just `np.random.seed(42)` if numpy already imported) after the imports section:
- `analysis/07_organoid_single_cell/02_qc.py`
- `analysis/07_organoid_single_cell/06_cell_cycle.py`
- `analysis/07_organoid_single_cell/07_single_cell_pathways.py`
- `analysis/07_organoid_single_cell/08_inhibitor_mechanism_exploration.py`

- [ ] **Step 6: Add seeds to section 07 preprocess.py**

In `analysis/07_organoid_single_cell/03_preprocess.py`:
- Add `np.random.seed(42)` after imports
- Line 236: change `sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)` to add `random_state=42`
- Line 338: change `sc.tl.umap(adata)` to `sc.tl.umap(adata, random_state=42)`
- Line 359: change `sc.tl.leiden(adata, resolution=resolution)` to `sc.tl.leiden(adata, resolution=resolution, random_state=42)`

- [ ] **Step 7: Commit**

```bash
git add analysis/03_rat_wes/01_parse_vep.py \
        analysis/03_rat_wes/02_cosmic_signatures.py \
        analysis/04_human_scrnaseq/16_run_cellphonedb.py \
        analysis/02_rat_snrnaseq/01_qc_filter.R \
        analysis/02_rat_snrnaseq/02_normalize_integrate.R \
        analysis/07_organoid_single_cell/02_qc.py \
        analysis/07_organoid_single_cell/03_preprocess.py \
        analysis/07_organoid_single_cell/06_cell_cycle.py \
        analysis/07_organoid_single_cell/07_single_cell_pathways.py \
        analysis/07_organoid_single_cell/08_inhibitor_mechanism_exploration.py
git commit -m "fix: add missing random seeds for reproducibility (W14,W15,W21,W28)"
```

---

## Task 3: Substantive Code Fixes (Batch 3)

Three changes that affect code logic.

**Files:**
- Modify: `analysis/04_human_scrnaseq/15_multicelltype_pathway.R:237-243`
- Modify: `analysis/01_human_bulk_rnaseq/01_preprocess.R:184-185,393-394,413-420`

### 3a. H1 — Unify curated pathway sets

- [ ] **Step 1: Replace FDR curated list with visualization list**

In `analysis/04_human_scrnaseq/15_multicelltype_pathway.R`, replace lines 237-243:

Old (FDR set):
```r
hallmark_curated <- c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_COMPLEMENT",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS",
  "LI_ESTROGENE_EARLY_E2_RESPONSE_UP", "LI_ESTROGENE_LATE_E2_RESPONSE_UP"
```

New (matches visualization set):
```r
hallmark_curated <- c(
  "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS",
  "LI_ESTROGENE_EARLY_E2_RESPONSE_UP", "LI_ESTROGENE_LATE_E2_RESPONSE_UP"
```

- [ ] **Step 2: Commit H1 fix**

```bash
git add analysis/04_human_scrnaseq/15_multicelltype_pathway.R
git commit -m "fix: unify FDR and visualization curated pathway sets (H1)"
```

### 3b. W24 — Duplicate gene removal by highest expression

- [ ] **Step 3: Fix TPM mode duplicate handling**

In `analysis/01_human_bulk_rnaseq/01_preprocess.R`, replace lines 184-185:

Old:
```r
  tpm_matrix_nodup <- tpm_matrix %>%
    dplyr::filter(!duplicated(GeneSymb))
```

New:
```r
  # Resolve duplicate gene symbols by keeping highest mean expression
  dup_genes <- unique(tpm_matrix$GeneSymb[duplicated(tpm_matrix$GeneSymb)])
  n_dups <- length(dup_genes)
  if (n_dups > 0) {
    cat("  Resolving", n_dups, "duplicated gene symbols (keeping highest mean expression)\n")
    keep <- rep(TRUE, nrow(tpm_matrix))
    expr_cols <- setdiff(colnames(tpm_matrix), "GeneSymb")
    for (g in dup_genes) {
      idx <- which(tpm_matrix$GeneSymb == g)
      means <- rowMeans(tpm_matrix[idx, expr_cols, drop = FALSE], na.rm = TRUE)
      keep[idx] <- FALSE
      keep[idx[which.max(means)]] <- TRUE
    }
    tpm_matrix_nodup <- tpm_matrix[keep, ]
  } else {
    tpm_matrix_nodup <- tpm_matrix
  }
```

- [ ] **Step 4: Fix count mode duplicate handling**

In `analysis/01_human_bulk_rnaseq/01_preprocess.R`, replace lines 393-394 with the same pattern adapted for the VST matrix:

Old:
```r
  vst_matrix_nodup <- vst_matrix %>%
    dplyr::filter(!duplicated(GeneSymb))
```

New:
```r
  # Resolve duplicate gene symbols by keeping highest mean expression
  dup_genes <- unique(vst_matrix$GeneSymb[duplicated(vst_matrix$GeneSymb)])
  n_dups <- length(dup_genes)
  if (n_dups > 0) {
    cat("  Resolving", n_dups, "duplicated gene symbols (keeping highest mean expression)\n")
    keep <- rep(TRUE, nrow(vst_matrix))
    expr_cols <- setdiff(colnames(vst_matrix), "GeneSymb")
    for (g in dup_genes) {
      idx <- which(vst_matrix$GeneSymb == g)
      means <- rowMeans(vst_matrix[idx, expr_cols, drop = FALSE], na.rm = TRUE)
      keep[idx] <- FALSE
      keep[idx[which.max(means)]] <- TRUE
    }
    vst_matrix_nodup <- vst_matrix[keep, ]
  } else {
    vst_matrix_nodup <- vst_matrix
  }
```

### 3c. M18 — Add annotation to count-mode heatmap

- [ ] **Step 5: Add annotation_col to distance heatmap**

In `analysis/01_human_bulk_rnaseq/01_preprocess.R`, around lines 415-417, add `annotation_col` to the `pheatmap()` call. First, build the annotation data frame from the sample annotation (same pattern as TPM mode around line 242):

```r
  # Build annotation for heatmap
  annot_for_heatmap <- data.frame(
    AgeRange = sample_annot$AgeRange,
    Group = sample_annot$Group,
    row.names = sample_annot$SampleName
  )

  pheatmap::pheatmap(
    sample_dist_matrix,
    main = "Sample Distance Matrix",
    annotation_col = annot_for_heatmap,
    ...  # keep existing params
  )
```

Check what variables are available in scope at that point — the annotation data frame may need to be constructed from `colData(dds)` or `sample_annot`.

- [ ] **Step 6: Commit W24 + M18 fixes**

```bash
git add analysis/01_human_bulk_rnaseq/01_preprocess.R
git commit -m "fix: duplicate gene handling by highest expression, annotate count-mode heatmap (W24, M18)"
```

---

## Task 4: VEP Cache Update + WES Re-run (Batch 4)

**Files:**
- Modify: `analysis/03_rat_wes/00_run_vep.sbatch:21,53`

- [ ] **Step 1: Download VEP cache v104**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/
curl -O https://ftp.ensembl.org/pub/release-104/variation/vep/rattus_norvegicus_vep_104_Rnor_6.0.tar.gz
tar xzf rattus_norvegicus_vep_104_Rnor_6.0.tar.gz
```

Expected: creates `rattus_norvegicus/104_Rnor_6.0/` alongside existing `rattus_norvegicus/95_Rnor_6.0/`.

- [ ] **Step 2: Update VEP sbatch**

In `analysis/03_rat_wes/00_run_vep.sbatch`:
- Line 53: change `--cache_version 95` to `--cache_version 104`
- Verify line 21 (`VEP_CACHE`) still points to the parent dir containing `rattus_norvegicus/`

- [ ] **Step 3: Back up current outputs**

```bash
cp -r analysis/03_rat_wes/outputs analysis/03_rat_wes/outputs_v95_backup
```

- [ ] **Step 4: Submit VEP re-annotation job**

Submit `00_run_vep.sbatch` via sbatch. Wait for completion.

- [ ] **Step 5: Submit analysis pipeline**

Submit `run_analysis.sbatch` (depends on VEP completion). This runs: parse VEP → COSMIC signatures → oncoplot.

- [ ] **Step 6: Spot-check results**

Compare `outputs/cosmic_signatures.csv` (new) vs `outputs_v95_backup/cosmic_signatures.csv` (old):
- Same dominant signatures per sample?
- Same total mutation counts?

Compare `outputs/oncoplot_data.csv` vs backup:
- Same genes present?
- Any consequence type changes?

If concordant: proceed. If materially different: stop and flag for user review.

- [ ] **Step 7: Commit**

```bash
git add analysis/03_rat_wes/00_run_vep.sbatch
git commit -m "fix: update VEP cache v95→v104 for Rnor_6.0 annotations (W6)"
```

---

## Task 5: MICA Integration (Batch 5)

Backfill the stub scripts in `analysis/01_human_bulk_rnaseq/mica/` with working code from `Jian_MICA/`.

**Files:**
- Modify: `analysis/01_human_bulk_rnaseq/mica/01_prep_data.R` (reconcile, not overwrite)
- Replace: `analysis/01_human_bulk_rnaseq/mica/03_run_mica.R`
- Replace: `analysis/01_human_bulk_rnaseq/mica/04_fig4_estrogen.R`
- Replace: `analysis/01_human_bulk_rnaseq/mica/05_fig5_immune.R`
- Create: `analysis/01_human_bulk_rnaseq/mica/00_setup_env.sh`
- Create: `analysis/01_human_bulk_rnaseq/mica/00_unzip_data.sh`
- Replace: `analysis/01_human_bulk_rnaseq/mica/02_run_gsva.R` (existing stub)
- Create: `analysis/01_human_bulk_rnaseq/mica/06_extract_claims.py`
- Create: `analysis/01_human_bulk_rnaseq/mica/07_statistical_audit.R`
- Create: `analysis/01_human_bulk_rnaseq/mica/08_compile_report.py`
- Modify: `analysis/01_human_bulk_rnaseq/mica/run_all.sbatch`

**Source directory:** `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/Jian_MICA/`

- [ ] **Step 1: Copy data files**

```bash
cp -r /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/Jian_MICA/data/* \
      analysis/01_human_bulk_rnaseq/mica/data/input/
```

Verify: `ls analysis/01_human_bulk_rnaseq/mica/data/input/` should show `clinical/`, `database/`, `gmt/`, `ilc_clean_data.RData`, `raw_data/`.

- [ ] **Step 2: Read and understand Jian_MICA scripts**

Read all scripts in `Jian_MICA/scripts/` (00-06) to understand:
- What each script does
- What paths they reference
- What environment they expect
- Whether they have `set.seed()`

Also read `Jian_MICA/required_data_files.md` for expected inputs.

- [ ] **Step 3: Copy and adapt setup scripts**

Copy `Jian_MICA/scripts/00_setup_env.sh` and `00_unzip_data.sh` to `mica/`. Update paths to reference `mica/data/input/` instead of `Jian_MICA/data/`.

- [ ] **Step 4: Reconcile 01_prep_data.R**

Read existing `mica/01_prep_data.R` (42 lines with `create_age_groups()` logic). Read `Jian_MICA/scripts/01_estrogene_gsva.R`. Merge:
- Keep the existing age group definitions if they match
- Add Jian's GSVA implementation
- Update all file paths to use `mica/data/input/` and `mica/outputs/`
- Ensure `set.seed(12345)` is present

- [ ] **Step 5: Create 02_run_gsva.R if needed**

If `Jian_MICA/scripts/01_estrogene_gsva.R` contains both data prep AND GSVA, the GSVA portion may need its own script `02_run_gsva.R`. If it's all in one script, merge into `01_prep_data.R` and skip this step.

- [ ] **Step 6: Replace 03_run_mica.R**

Copy `Jian_MICA/scripts/02_mica_analysis.R` to `mica/03_run_mica.R`. Update:
- Input paths: `mica/outputs/` for GSVA results
- Output paths: `mica/outputs/` for MICA results
- Add `set.seed(12345)` if missing

- [ ] **Step 7: Replace figure generation scripts**

Read `Jian_MICA/scripts/03_figure_generation.R`. Split into:
- `mica/04_fig4_estrogen.R` — Figure 4C (HSD17B heatmap) and 4H (EstroGene boxplots)
- `mica/05_fig5_immune.R` — Figure 5C (inflammatory heatmap) and 5E (immune boxplots)

Merge helper functions from existing stubs (`heatmap_median()`, `path_boxplot()`, `cell_boxplot()`) with Jian's implementations. Use Jian's version if they differ — it's the one that produced the published figures.

Update output paths to `mica/figures/`.

- [ ] **Step 8: Adapt remaining scripts**

Copy and adapt:
- `Jian_MICA/scripts/04_extract_claims.py` → `mica/06_extract_claims.py`
- `Jian_MICA/scripts/05_statistical_audit.R` → `mica/07_statistical_audit.R`
- `Jian_MICA/scripts/06_compile_report.py` → `mica/08_compile_report.py`

Update all paths. Add seeds where missing.

- [ ] **Step 9: Create unified sbatch**

Replace `mica/run_all.sbatch` with a unified script that:
- Uses `erp_brca_aging` conda environment
- Includes `#SBATCH --mail-type=FAIL` and `#SBATCH --mail-user=alc376@pitt.edu`
- Runs scripts 01-08 sequentially with `set -eo pipefail`
- Writes completion marker to `.pipeline_markers/01_mica.complete`

- [ ] **Step 10: Verify execution**

Submit `mica/run_all.sbatch`. Compare outputs against `Jian_MICA/results/` and `Jian_MICA/figures/`:
- Do `fig4c.png`, `fig4h.png`, `fig5c.png`, `fig5e.png` match?
- Do MICA scheme1/scheme2 CSVs match?

- [ ] **Step 11: Commit**

```bash
git add analysis/01_human_bulk_rnaseq/mica/
git commit -m "feat: integrate MICA concordance analysis from Jian_MICA into monorepo (C3,W10,W25)"
```

---

## Task 6: Methods Documentation Corrections (Batch 6)

Depends on Tasks 1-5 being complete. Update all methods files to match the corrected code.

**Files:**
- Modify: `docs/methods/01_human_bulk_rnaseq_methods.md`
- Modify: `docs/methods/02_rat_snrnaseq_methods.md`
- Modify: `docs/methods/03_rat_wes_methods.md`
- Modify: `docs/methods/04_human_scrnaseq_methods.md`
- Modify: `docs/methods/05_rat_bulk_rnaseq_methods.md`
- Modify: `docs/methods/06_spatial_biopsies_methods.md`
- Modify: `docs/methods/07_organoid_scrnaseq_methods.md`
- Modify: `docs/methods/00_main_text_methods.md`

- [ ] **Step 1: Update section 01 methods**

In `docs/methods/01_human_bulk_rnaseq_methods.md`:
- Line 83-84: Remove "Log2(TPM) expression matrix" from PROGENy input. State it uses `vst_normalized_matrix.rds` (H6)
- Section 1 (Age Groups): Add numeric boundaries for primary analysis age groups — check the annotation file to determine actual values (M20/W22)
- Section 4.1: Change "Planned" to describe the executed MICA analysis with actual parameters and results
- Line 44-45: Update duplicate gene handling from "first occurrence" to "highest mean expression across samples"
- Update main text summary at bottom to match

- [ ] **Step 2: Update section 02 methods**

In `docs/methods/02_rat_snrnaseq_methods.md`:
- Section 6 (DE): Change "Cell types with fewer than 10 cells in either group were excluded" to "Cell types with fewer than 10 cells per sample were excluded; groups required at least 2 samples" (H2)
- Section 6: Add "Gene filtering required >= 10 counts in >= 2 samples before DESeq2" (M13)
- Section 5/6: State explicitly "Differential expression and differential abundance analyses used the marker-scoring annotation from the primary pipeline (Louvain resolution 0.4), not scType annotations" (H3/W11)
- Section 10 (SLURM): Change `erp_brca_aging` to `erp_snrnaseq` (M1)
- Section 7: Add note explaining why propeller's p-values are re-corrected with `p.adjust()` (M12)
- Section 2.2: Document that DoubletFinder PCA computes 50 PCs (Seurat default), with dims 1:30 used for paramSweep and doubletFinder (M4)

- [ ] **Step 3: Update section 03 methods**

In `docs/methods/03_rat_wes_methods.md`:
- Section 2: Change `--cache_version 95` to `--cache_version 104`
- Section 6 (Pipeline): Change conda env from `erp_brca_aging` to `aging_wes` (M2)
- Section 8 (Reproducibility): Add note: "Results were regenerated with VEP v114.2 using the Ensembl release 104 cache (latest available for Rnor_6.0). The original analysis used VEP v95 with cache v95."

- [ ] **Step 4: Update section 04 methods**

In `docs/methods/04_human_scrnaseq_methods.md`:
- Section 23 (Software): Change R version from `4.4.1` to `4.3.3` (H5 — verified)
- Section 18: Document that script 15 uses `min_patients_per_group = 2` (vs script 14's `MIN_PATIENTS_PER_GROUP = 3`) (M10)
- Section 22: Change "top50" mode description to note it actually selects 40 pairs (M11)
- Section 16: Clarify that script 13 generates a bar chart of communication-related DEGs, not an L-R interaction analysis (M17)
- Section 18: Update curated pathway set description to match the unified set (after H1 fix — now includes IL2_STAT5, IFN_ALPHA, EMT; excludes COMPLEMENT)
- Add note that script 09 GSEA uses Wilcoxon z-score for gene ranking (not log2FC * -log10(p)) and queries 3 libraries: MSigDB_Hallmark_2020, KEGG_2021_Human, Reactome_2022 (M14-15)

- [ ] **Step 5: Update section 05 methods**

In `docs/methods/05_rat_bulk_rnaseq_methods.md`:
- Section 7.2: Fix validation matching description. HTSeq count comparison uses a hardcoded sample name lookup table (not gsub). PAM50 comparison uses `gsub("-Tumor", "", ...)` for fuzzy matching. (M19)

- [ ] **Step 6: Update section 06 methods**

In `docs/methods/06_spatial_biopsies_methods.md`:
- Section 4: Change "Results from enrichment and prerank analyses were combined" to "Separate summary heatmaps were generated for enrichment and prerank results" (M16)
- Add CITEgeist version or tool citation (V2)
- Add full Carleton et al. citation with year, journal, DOI (V3)

- [ ] **Step 7: Update section 07 methods**

In `docs/methods/07_organoid_scrnaseq_methods.md`:
- Section 2: Add "Cell Ranger v9.0.1" (M5)
- Section 2: Add reference build details: "Custom reference built from Ensembl release 109 primary assembly FASTA and GENCODE v44 GTF, with biotype filtering and PAR_Y gene exclusion" (M6)
- Section 7: Add note distinguishing script 05 GSEA (pseudobulk-based, log2FC * -log10(p) ranking, Hallmarks only) from script 09 GSEA (cell-level Wilcoxon z-score ranking, 3 libraries: Hallmark, KEGG, Reactome) (M14-15)
- Section 6: Add acknowledgment: "Note that with one biological replicate per treatment condition, the pseudobulk DE analysis has limited statistical power; results should be interpreted as exploratory" (V6)

- [ ] **Step 8: Update main text methods**

In `docs/methods/00_main_text_methods.md`:
- Section 01 paragraph: Update to reflect MICA concordance analysis as executed (not planned)
- Section 03 paragraph: Update VEP cache to v104

- [ ] **Step 9: Commit methods updates**

```bash
git add docs/methods/
git commit -m "docs: update all methods sections to match corrected code"
```

---

## Task 7: Update Code Review Issues Tracker

**Files:**
- Modify: `docs/code_review/code_review_issues.md`

- [ ] **Step 1: Read current issues file**

Read `docs/code_review/code_review_issues.md` to get current state.

- [ ] **Step 2: Remove W18**

Delete the W18 row entirely (wrong project — ERT Visium duplicate sample, confirmed not in this repo's data).

- [ ] **Step 3: Mark resolved issues with strikethrough**

Apply `~~strikethrough~~` + `**RESOLVED**` pattern (matching existing resolved entries) to:

**Critical:** C3 (MICA stubs → integrated)

**Warnings:** W3, W6, W10, W11, W13, W14, W15, W20, W21, W24, W25, W26, W27, W28, W29, W35, W36

For each, add a brief resolution note matching the pattern of existing resolved entries (e.g., `**RESOLVED** — comment fixed to match code`).

- [ ] **Step 4: Add new resolved entry for H1**

Add a new warning entry (W37 or next available number):
```
| W37 | 04_human_scrnaseq | ~~FDR and visualization curated pathway sets diverged~~ — `15_multicelltype_pathway.R` used different pathway lists for FDR correction (lines 237-243) vs visualization (lines 393-401). FDR set included COMPLEMENT but excluded IL2_STAT5, IFN_ALPHA, EMT. | **RESOLVED** — FDR set unified with visualization set (15 Hallmark + 12 BioCarta) |
```

- [ ] **Step 5: Commit**

```bash
git add docs/code_review/code_review_issues.md
git commit -m "docs: update code review issues tracker — resolve 19 issues, remove W18"
```

---

## Execution Summary

| Task | Batch | Parallelizable | HPC Required | Dependencies |
|------|-------|---------------|-------------|--------------|
| 1 (Quick fixes) | 1 | Yes (with 2,3) | No | None |
| 2 (Seeds) | 2 | Yes (with 1,3) | No | None |
| 3 (Substantive) | 3 | Yes (with 1,2) | No | None |
| 4 (VEP + WES) | 4 | Yes (with 1-3) | Yes (sbatch) | Cache download |
| 5 (MICA) | 5 | Yes (with 1-4) | Yes (sbatch to verify) | Data copy |
| 6 (Methods docs) | 6 | No | No | Tasks 1-5 |
| 7 (Issues tracker) | 6 | No | No | Tasks 1-5 |
