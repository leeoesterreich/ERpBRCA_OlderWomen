# Supplementary Methods

## Table of Contents

1. [Human Bulk RNA-seq Analysis](#1-human-bulk-rna-seq-analysis)
2. [Rat snRNA-seq Analysis](#2-rat-snrna-seq-analysis)
3. [Rat Whole Exome Sequencing](#3-rat-whole-exome-sequencing)
4. [Human Single-Cell RNA-seq Analysis](#4-human-single-cell-rna-seq-analysis)
5. [Rat Bulk RNA-seq Analysis](#5-rat-bulk-rna-seq-analysis)
6. [Issues Requiring Resolution](#6-issues-requiring-resolution)

---

For detailed per-subsection methods with inline issue flags, see individual files:
- `01_human_bulk_rnaseq_methods.md`
- `02_rat_snrnaseq_methods.md`
- `03_rat_wes_methods.md`
- `04_human_scrnaseq_methods.md`
- `05_rat_bulk_rnaseq_methods.md`

---

## 6. Issues Requiring Resolution

The following issues were identified during code review and methods extraction. Items marked **CRITICAL** indicate code behavior that produces incorrect or misleading results and must be fixed before the methods can be accurately reported. Items marked **WARNING** indicate concerns that should be reviewed and documented.

### Critical Issues

| # | Subsection | Issue | Impact |
|---|-----------|-------|--------|
| C1 | 04_human_scrnaseq | ~~DEG column schema mismatch~~ | **RESOLVED** — scripts 12-13 updated to use DESeq2 column names (`log2FoldChange`, `padj`) |
| C2 | 04_human_scrnaseq | ~~Multi-cell-type scripts load macrophage-only data~~ | **RESOLVED** — scripts 15-16 now load `seurat_annotated.rds` (full multi-cell-type object) |
| C3 | 01_human_bulk_rnaseq | MICA sub-pipeline (5 scripts) contains only stub implementations with TODO comments | No MICA results were generated; any reported MICA findings are not from this pipeline |
| C4 | 01_human_bulk_rnaseq | ~~Correlation analysis uses unfiltered gene set~~ | **RESOLVED** — `04_correlations.R` now filters to protein-coding genes (consistent with `03_run_progeny.R`) |

### Warnings

| # | Subsection | Issue |
|---|-----------|-------|
| W1 | 02_rat_snrnaseq | ~~Cell-level DE pseudoreplication~~ | **RESOLVED** — converted to pseudobulk DESeq2 with per-sample aggregation |
| W2 | 05_rat_bulk_rnaseq | ~~CPM mislabeled as TPM~~ | **RESOLVED** — renamed to `normalized_cpm.csv`, updated all references |
| W3 | 04_human_scrnaseq | Script 15 comment says Wilcoxon but code runs Welch t-test |
| W4 | 04_human_scrnaseq | Dense matrix conversion in GSVA (~50 GB) risks OOM on standard nodes |
| W5 | 04_human_scrnaseq | ~~CellPhoneDB p-values unadjusted~~ | **RESOLVED** — BH-FDR correction now applied before filtering/display |
| W6 | 03_rat_wes | VEP cache v95 paired with VEP software v114.2 — major version mismatch |
| W7 | 03_rat_wes | SigProfilerAssignment and pybiomart not in environment.yml; separate untracked env used |
| W8 | 03_rat_wes | BioMart live API dependency (no archive host pinned) introduces non-reproducibility |
| W9 | 01_human_bulk_rnaseq | ~~PROGENy reads raw TPM independently~~ | **RESOLVED** — PROGENy now uses same VST matrix as GSVA |
| W10 | 01_human_bulk_rnaseq | Two different conda environments referenced across sbatch scripts |
| W11 | 02_rat_snrnaseq | Annotation ambiguity: primary vs scType pipeline outputs diverge; unclear which feeds DE/DA |
| W12 | 02_rat_snrnaseq | Fixed 5% doublet rate not calibrated to loading density |
| W13 | 03_rat_wes | Misleading comment in oncoplot script swaps Young/Old label descriptions |
| W14 | 03_rat_wes, 04_human_scrnaseq | Missing random seeds in Python scripts (project requirement violation) |
| W15 | 04_human_scrnaseq | Missing random seed in CellPhoneDB v5 script (v4 uses seed 42) |
| W16 | 05_rat_bulk_rnaseq | BioMart version not pinned in DESeq2 and PAM50 scripts |
| W17 | 05_rat_bulk_rnaseq | STAR genome index provenance undocumented |
| W18 | 01_human_bulk_rnaseq | Duplicate sample 135938ERPRnoNAC not explicitly excluded (per CLAUDE.md requirement) |
| W19 | All | Most package versions unpinned in environment.yml; no runtime sessionInfo() captured |
