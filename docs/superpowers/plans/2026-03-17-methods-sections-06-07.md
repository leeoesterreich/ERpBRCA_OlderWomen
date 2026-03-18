# Methods Documentation (Sections 06 & 07 + Cleanup) Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write clean manuscript-ready supplementary methods for sections 06 (spatial biopsies) and 07 (organoid single-cell), clean internal review notes from existing methods (01-05), and consolidate issues into a separate document.

**Architecture:** This is a documentation-only plan — no code changes, no tests. Each task produces one or two markdown files following the existing prose style in `docs/methods/`. Tasks are independent and can be parallelized. The existing methods files use numbered sections with parameterized prose, embedded tool versions, and reproducibility notes.

**Tech Stack:** Markdown, git

**Spec:** `docs/superpowers/specs/2026-03-17-methods-sections-06-07-design.md`

---

## File Map

| Action | File | Responsibility |
|--------|------|---------------|
| CREATE | `docs/methods/06_spatial_biopsies_methods.md` | New supplementary methods for spatial biopsies |
| CREATE | `docs/methods/07_organoid_scrnaseq_methods.md` | New supplementary methods for organoid scRNA-seq |
| CREATE | `docs/code_review/code_review_issues.md` | Consolidated issues extracted from all methods docs |
| MODIFY | `docs/methods/00_main_text_methods.md` | Add paragraphs for sections 06 and 07 |
| MODIFY | `docs/methods/00b_supplementary_methods.md` | Clean to TOC-only, add sections 06-07 |
| MODIFY | `docs/methods/01_human_bulk_rnaseq_methods.md` | Strip warnings/issues, keep methods |
| MODIFY | `docs/methods/02_rat_snrnaseq_methods.md` | Strip warnings/issues, keep methods |
| MODIFY | `docs/methods/03_rat_wes_methods.md` | Strip warnings/issues, keep methods |
| MODIFY | `docs/methods/04_human_scrnaseq_methods.md` | Strip warnings/issues, keep methods |
| MODIFY | `docs/methods/05_rat_bulk_rnaseq_methods.md` | Strip warnings/issues, keep methods |

---

### Task 1: Write section 06 supplementary methods

**Files:**
- Create: `docs/methods/06_spatial_biopsies_methods.md`

**Context for writer:** Read the spec section for `06_spatial_biopsies_methods.md` (lines 10-34 of the design spec). Also read `docs/methods/05_rat_bulk_rnaseq_methods.md` as a style reference — use numbered sections, flowing prose with embedded parameters, not bullet lists. Cross-reference against analysis scripts `analysis/06_spatial_biopsies/01_immune_secretion.py` and `analysis/06_spatial_biopsies/02_immune_pathways.py` for accuracy.

- [ ] **Step 1: Read style reference and analysis scripts**

Read `docs/methods/05_rat_bulk_rnaseq_methods.md` for prose style. Read `analysis/06_spatial_biopsies/01_immune_secretion.py` and `analysis/06_spatial_biopsies/02_immune_pathways.py` for exact parameters. Read `analysis/06_spatial_biopsies/figure_config.py` for figure settings.

- [ ] **Step 2: Write `06_spatial_biopsies_methods.md`**

Create `docs/methods/06_spatial_biopsies_methods.md` with these sections:
1. Data Source & Deconvolution
2. Spatial Cytokine Visualization
3. Macrophage Pathway Enrichment (enrichr + GSEA prerank)
4. Summary Visualization
5. Software & Reproducibility

Follow spec lines 14-34. Use numbered sections and flowing prose. Include exact parameter values from the scripts.

- [ ] **Step 3: Verify against spec**

Confirm all items from spec lines 14-34 are covered. Check that the dotplot description correctly notes the interleukin-only filtering (excluding TGF terms at visualization stage).

- [ ] **Step 4: Commit**

```bash
git add docs/methods/06_spatial_biopsies_methods.md
git commit -m "docs: add supplementary methods for spatial biopsies (section 06)"
```

---

### Task 2: Write section 07 supplementary methods

**Files:**
- Create: `docs/methods/07_organoid_scrnaseq_methods.md`

**Context for writer:** Read the spec section for `07_organoid_scrnaseq_methods.md` (lines 38-122 of the design spec). Use `docs/methods/05_rat_bulk_rnaseq_methods.md` as style reference. Cross-reference against `analysis/07_organoid_single_cell/_config.py` for constants (seeds, thresholds, paths), and scripts `02_qc.py` through `11_proliferation_analysis.py` for exact parameters.

- [ ] **Step 1: Read style reference, _config.py, and analysis scripts**

Read `docs/methods/05_rat_bulk_rnaseq_methods.md` for prose style. Read `analysis/07_organoid_single_cell/_config.py` for shared constants. Read key scripts: `02_qc.py`, `03_preprocess.py`, `04_pseudobulk_de.py`, `05_pathways.py`, `06_cell_cycle.py`, `07_single_cell_pathways.py`, `08_inhibitor_mechanism_exploration.py`, `09_pathway_enrichment_analysis.py`, `10_heterogeneity_analysis.py`, `11_proliferation_analysis.py`.

- [ ] **Step 2: Write `07_organoid_scrnaseq_methods.md`**

Create `docs/methods/07_organoid_scrnaseq_methods.md` with these sections:
1. Experimental Design
2. Cell Ranger & Reference (include filter-probes=false rationale)
3. Quality Control & Doublet Detection
4. Preprocessing & Dimensionality Reduction
5. Cell Cycle Scoring
6. Pseudobulk Differential Expression
7. Pathway & Transcription Factor Analysis (GSEA prerank + Enrichr + PROGENy + DoRothEA)
8. Cell Cycle Distribution Analysis (chi-square)
9. Single-Cell Pathway Scoring
10. HSD17B7 Inhibitor Mechanism Analysis
11. Heterogeneity Analysis (including estrogen response continuum)
12. Proliferation & Quiescence Classification
13. Treatment Label Convention
14. Software & Reproducibility

Follow spec lines 42-122. Use numbered sections and flowing prose. Note seed=42 exception and the distinction between pseudobulk DE and cell-level Wilcoxon DE.

- [ ] **Step 3: Verify against spec**

Confirm all items from spec lines 42-122 are covered. Verify: G2M gene count is 54 (not 47), seed documentation notes which scripts are seeded, proliferation lists all 7 treatment pairs, Enrichr lists all 5 libraries.

- [ ] **Step 4: Commit**

```bash
git add docs/methods/07_organoid_scrnaseq_methods.md
git commit -m "docs: add supplementary methods for organoid scRNA-seq (section 07)"
```

---

### Task 3: Create consolidated issues document

**Files:**
- Create: `docs/code_review/code_review_issues.md`

**Context for writer:** Extract all issues from two sources: (1) `docs/methods/00b_supplementary_methods.md` lines 23-58 (the existing issue tracker with C1-C4 and W1-W19), and (2) any new issues found in sections 06-07 (continue numbering from W20+, C5+). Read the existing methods docs (01-05) for any inline warnings (`> **[WARNING]**`, `> **[CRITICAL]**`, `> **[RESOLVED]**`) that should also be captured.

- [ ] **Step 1: Read existing issues from `00b_supplementary_methods.md`**

Read `docs/methods/00b_supplementary_methods.md` to get C1-C4 and W1-W19.

- [ ] **Step 2: Scan existing methods files for inline issues**

Read `docs/methods/01_human_bulk_rnaseq_methods.md` through `05_rat_bulk_rnaseq_methods.md`. Search for blockquote patterns (`> **[WARNING]**`, `> **[CRITICAL]**`, `> **[RESOLVED]**`) and any inline `⚠️` markers. Catalog every issue not already in C1-C4/W1-W19.

- [ ] **Step 3: Identify new issues from sections 06-07**

Read `analysis/06_spatial_biopsies/02_immune_pathways.py` and `analysis/07_organoid_single_cell/` scripts. Known new issues to include:
- W20: gseapy `organism='Human'` parameter in section 06 (deprecated/ignored in recent versions)
- W21: Section 07 scripts 02, 03, 06, 07, 08 do not set explicit random seeds (project convention is seed=42 for this section)
- Add any others found during scanning.

- [ ] **Step 4: Write `docs/code_review/code_review_issues.md`**

Create the file organized by section (01-07). Each issue gets: ID, section, severity (CRITICAL/WARNING), description, resolution status. Preserve all existing issue text and resolution notes. Use a table format matching the existing `00b` style.

- [ ] **Step 5: Commit**

```bash
git add docs/code_review/code_review_issues.md
git commit -m "docs: create consolidated code review issues document"
```

---

### Task 4: Clean existing methods files (01-05)

**Files:**
- Modify: `docs/methods/01_human_bulk_rnaseq_methods.md`
- Modify: `docs/methods/02_rat_snrnaseq_methods.md`
- Modify: `docs/methods/03_rat_wes_methods.md`
- Modify: `docs/methods/04_human_scrnaseq_methods.md`
- Modify: `docs/methods/05_rat_bulk_rnaseq_methods.md`

**Context for writer:** For each file, remove all internal review material while keeping the analytical methods content. The consolidated issues document (Task 3) now holds all issue tracking.

**What to remove:**
- Entire blockquote blocks starting with `> **[WARNING]**`, `> **[CRITICAL]**`, `> **[RESOLVED]**` (remove the full blockquote, not just the tag)
- Inline `⚠️` markers and their associated commentary sentences
- Sections about planned-but-unimplemented features (MICA stubs in section 01, MCPcounter, immune deconvolution that was never implemented)
- Bug workaround notes (Seurat 5.x patchwork, applymap deprecation, etc.)
- "Critical Issues" or "Known Issues" sections at the top of files

**What to keep:**
- All parameter values, statistical methods, thresholds, tool versions
- Analytical workflow descriptions
- Reproducibility notes (seeds, determinism)
- Section structure and numbered headers

- [ ] **Step 1: Clean `01_human_bulk_rnaseq_methods.md`**

Read file. Remove all blockquote warnings, inline issue commentary, MICA stub documentation, MCPcounter/immune deconvolution plans. Keep all methods content.

- [ ] **Step 2: Clean `02_rat_snrnaseq_methods.md`**

Read file. Remove all blockquote warnings (doublet rate calibration, scType issues, Harmony seed, etc.). Keep methods content.

- [ ] **Step 3: Clean `03_rat_wes_methods.md`**

Read file. Remove all blockquote warnings (VEP cache mismatch, applymap deprecation, missing seeds, BioMart API). Keep methods content.

- [ ] **Step 4: Clean `04_human_scrnaseq_methods.md`**

Read file. This is the largest file (281 lines). Remove the "Critical Issues" section at top, all inline warnings (Wilcoxon/Welch mismatch, OOM risk, CellPhoneDB seeds, hardcoded paths). Keep methods content.

- [ ] **Step 5: Clean `05_rat_bulk_rnaseq_methods.md`**

Read file. Remove blockquote about CPM/TPM rename (line 31 area), any other warnings. Keep methods content.

- [ ] **Step 6: Verify each file still reads as clean manuscript prose**

Skim each cleaned file. Ensure no orphaned `>` characters, no dangling references to removed content, no `⚠️` symbols remaining.

- [ ] **Step 7: Commit**

```bash
git add docs/methods/01_human_bulk_rnaseq_methods.md docs/methods/02_rat_snrnaseq_methods.md docs/methods/03_rat_wes_methods.md docs/methods/04_human_scrnaseq_methods.md docs/methods/05_rat_bulk_rnaseq_methods.md
git commit -m "docs: clean internal review notes from methods files (01-05)"
```

---

### Task 5: Update main text and supplementary TOC

**Files:**
- Modify: `docs/methods/00_main_text_methods.md`
- Modify: `docs/methods/00b_supplementary_methods.md`

**Context for writer:** Add two new paragraphs to the main text methods matching the existing dense one-paragraph-per-section style. Clean the supplementary TOC to remove the issue tracker and add sections 06-07.

- [ ] **Step 1: Read `00_main_text_methods.md` for style reference**

Read the file. Note the dense paragraph format with embedded parameters, citations, and tool versions.

- [ ] **Step 2: Add spatial biopsies paragraph**

After the "### Rat Bulk RNA-seq" section, add:

```markdown
### Spatial Biopsy Analysis
```

Write a dense paragraph covering: 6 Visium biopsies from HCC22-088 (Carleton et al.), CITEgeist deconvolution into 9 cell types (3 macrophage subtypes: CD163+, CD14+HLA-DR+, CD11c+), spatial cytokine mapping (IL4/IL10/IL13/TGFB1-3), macrophage pathway enrichment via gseapy (MSigDB Hallmark, KEGG, Reactome; top 200 genes, FDR<0.05) and GSEA prerank (1,000 permutations, seed 12345), interleukin pathway focus. Match spec lines 160-161.

- [ ] **Step 3: Add organoid scRNA-seq paragraph**

After the spatial biopsies section, add:

```markdown
### Organoid Single-Cell RNA-seq
```

Write a dense paragraph covering: PDO-296 (patient-derived from study cohort), 7 conditions (Vehicle, E1, E2, E1+fulvestrant, E2+fulvestrant, E1+HSD17B7i, E2+HSD17B7i), 10x Chromium Flex with custom GRCh38-2024-A reference (probe filtering disabled to retain HSD17B7), Scrublet QC (expected rate 0.06), CPM normalization + log1p, 3,000 HVGs, Leiden clustering (resolution 0.5), PyDESeq2 pseudobulk DE, GSEA prerank (Hallmark, seed 42), PROGENy/DoRothEA via decoupler, custom gene set scoring (8 sets, Mann-Whitney U, Cohen's d, BH-FDR), chi-square cell cycle distribution, GMM cycling classification, E1/E2-dominant signature analysis for inhibitor mechanism, heterogeneity (neighborhood mixing, estrogen continuum). Match spec lines 163-164.

- [ ] **Step 4: Clean `00b_supplementary_methods.md`**

Replace entire contents with a clean TOC pointing to all 7 sections:

```markdown
# Supplementary Methods

## Table of Contents

1. [Human Bulk RNA-seq Analysis](01_human_bulk_rnaseq_methods.md)
2. [Rat snRNA-seq Analysis](02_rat_snrnaseq_methods.md)
3. [Rat Whole Exome Sequencing](03_rat_wes_methods.md)
4. [Human Single-Cell RNA-seq Analysis](04_human_scrnaseq_methods.md)
5. [Rat Bulk RNA-seq Analysis](05_rat_bulk_rnaseq_methods.md)
6. [Spatial Biopsy Analysis](06_spatial_biopsies_methods.md)
7. [Organoid Single-Cell RNA-seq Analysis](07_organoid_scrnaseq_methods.md)

---

For code review issues and known limitations, see [`docs/code_review/code_review_issues.md`](../code_review/code_review_issues.md).
```

- [ ] **Step 5: Commit**

```bash
git add docs/methods/00_main_text_methods.md docs/methods/00b_supplementary_methods.md
git commit -m "docs: add main text methods for sections 06-07, clean supplementary TOC"
```

---

## Task Dependencies

Tasks 1-3 are fully independent and can be parallelized.

Task 4 depends on Task 3 (issues must be extracted before removing them from methods files).

Task 5 is independent of Tasks 1-4 (modifies different files).

```
Task 1 (section 06 methods)  ──┐
Task 2 (section 07 methods)  ──┼── can run in parallel
Task 3 (issues doc)          ──┤
Task 5 (main text + TOC)     ──┘
                                │
                                ▼
Task 4 (clean 01-05)         ── depends on Task 3
```
