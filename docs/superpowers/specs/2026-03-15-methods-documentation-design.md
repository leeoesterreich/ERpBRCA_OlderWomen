# Methods Documentation Design

**Date**: 2026-03-15
**Project**: ERpBRCA_OlderWomen (Carleton et al. 2024)
**Purpose**: Generate Nature-style Methods (main text + supplementary) from code

## Scope

Document the computational methods for all 5 analysis subsections by reading
every analysis script and extracting what the code actually does.

## Output Structure

```
docs/methods/
├── 00_main_text_methods.md           # Compressed Nature main-text (~500-800 words)
├── 00b_supplementary_methods.md      # Combined supplementary (all subsections)
├── 01_human_bulk_rnaseq_methods.md
├── 02_rat_snrnaseq_methods.md
├── 03_rat_wes_methods.md
├── 04_human_scrnaseq_methods.md
└── 05_rat_bulk_rnaseq_methods.md
```

## Per-subsection template

Each supplementary methods file follows this structure:

1. **Data acquisition** — GEO accession, sample counts, species, tissue type
2. **Preprocessing** — QC thresholds, filtering criteria, normalization method
3. **Core analysis** — Tools, exact parameters, statistical tests, MTC method
4. **Downstream analysis** — Pathway enrichment, cell-cell communication, etc.
5. **Visualization** — Plot types, figure panel mapping
6. **Software versions** — From environment.yml (ground truth)
7. **Reproducibility notes** — Seeds, checkpoints, deterministic ordering

## Issue flagging convention

Issues identified in code review are flagged inline using admonitions:

> **[CRITICAL]** Description of issue that makes this method suspect.
> The code should be revised before publication.

> **[WARNING]** Description of concern (e.g., pseudoreplication, naming mismatch).
> Review before finalizing methods text.

## Execution

5 parallel subagents, one per subsection. Each reads:
- All scripts in the subsection directory
- environment.yml for versions
- Writes the supplementary methods file
- Writes a compressed paragraph for main-text assembly

After all complete, assemble 00_main_text_methods.md and 00b_supplementary_methods.md.

## Known issues to flag (from Codex review)

### Critical
- 04: DEG column schema mismatch (DESeq2 vs Seurat style) breaks downstream
- 04: Multi-cell-type scripts fed macrophage-only data
- 01: Duplicate sample 135938ERPRnoNAC not excluded

### Warning
- 02: Cell-level DE (FindMarkers) instead of pseudobulk — pseudoreplication
- 05: "TPM" normalization lacks gene-length correction
- 04: Comment says Wilcoxon, code does Welch t-test
- 03/04: Missing random seeds in Python scripts
- 04: Dense matrix conversion OOM risk
- 04: CellPhoneDB p-values not adjusted at plot stage
