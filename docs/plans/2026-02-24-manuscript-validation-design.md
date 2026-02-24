# Manuscript Figure Validation Pipeline Design

**Date:** 2026-02-24
**Status:** Approved
**Author:** Claude (with user approval)

## Overview

Systematic validation that refactored codebase produces identical scientific conclusions to the published manuscript figures. Uses tiered comparison (numerical, visual, claims) with Azure OpenAI GPT-5 for figure analysis.

## Problem

The ERpBRCA_OlderWomen codebase has been refactored and corrected. Before publication/submission, we need to verify that all computational figures and claims in the manuscript are reproducible from the cleaned code.

## Scope

### In Scope (Computational Figures)

| Figure | Description | Analysis Folder |
|--------|-------------|-----------------|
| Fig. 2 | WES mutations, COSMIC signatures | `03_rat_wes` |
| Fig. 4 (partial) | HSD17B7 correlations with estrogen pathways | `01_human_bulk_rnaseq` |
| Fig. 6 (panels B-D) | scRNA-seq analysis | `04_human_scrnaseq` |
| Fig. 7 | MICA deconvolution + scRNA-seq TAMs | `01_human_bulk_rnaseq`, `04_human_scrnaseq` |
| EDF 2 | Rat multi-omic validation | `02_rat_snrnaseq`, `03_rat_wes`, `05_rat_bulk_rnaseq` |
| EDF 5B | MICA enzyme analysis | `01_human_bulk_rnaseq` |
| EDF 6 | HSD17B7 pathway correlations | `01_human_bulk_rnaseq` |
| EDF 7 | ESR1 age expression | `01_human_bulk_rnaseq` |
| EDF 9 | TME inflammation (if computational) | `01_human_bulk_rnaseq` |
| EDF 10 (partial) | Spatial transcriptomics | TBD |

### Out of Scope (Wet Lab Figures)

- Fig. 1: Rat tumor model, survival curves
- Fig. 3: Plasma E1/E2 measurements (Luminex/ELISA)
- Fig. 5: PDO inhibitor experiments
- Fig. 8: Macrophage IF experiments
- EDF 1: Rat estrogen/menstruation
- EDF 3-4: Tissue estrogen disposition
- EDF 8: PDO conversion assays
- Supplementary Tables S1-S3

## Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                    MANUSCRIPT VALIDATION PIPELINE                    │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────────────┐  │
│  │  EXTRACTION  │───▶│  COMPARISON  │───▶│  TIERED REPORT       │  │
│  └──────────────┘    └──────────────┘    └──────────────────────┘  │
│         │                   │                                        │
│         ▼                   ▼                                        │
│  ┌────────────┐      ┌────────────┐                                 │
│  │ PowerPoint │      │ Numerical  │  ← CSV/RDS programmatic compare │
│  │ → PNG      │      │ Tier       │                                 │
│  ├────────────┤      ├────────────┤                                 │
│  │ DOCX       │      │ Visual     │  ← Azure OpenAI GPT-5 vision    │
│  │ → Legends  │      │ Tier       │                                 │
│  │ → Claims   │      ├────────────┤                                 │
│  └────────────┘      │ Claim      │  ← Text match to data outputs   │
│                      │ Tier       │                                 │
│                      └────────────┘                                 │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

## Directory Structure

```
validation/
├── extracted/
│   ├── figures/
│   │   ├── main/           # PNGs from main PowerPoint
│   │   └── supplemental/   # PNGs from supplemental PowerPoint
│   ├── legends.json        # Structured figure legends
│   ├── claims.json         # Quantitative claims from manuscript text
│   └── figure_mapping.json # Maps figures to analysis outputs
├── comparisons/
│   ├── numerical/          # CSV diffs, correlation checks
│   ├── visual/             # Azure OpenAI analysis outputs
│   └── claims/             # Claim verification results
└── report/
    ├── validation_report.md
    └── figures/            # Side-by-side comparison images
```

## Input Files

| File | Location | Content |
|------|----------|---------|
| Main figures | `Manuscript/Aging_Main Figures_2-1-26.pptx` | 8 slides, Figs 1-8 |
| Supplemental | `Manuscript/Aging_Supplementary Material_2-1-26.pptx` | 19 slides, EDF 1-10, Supp Figs 1-6, Tables S1-S3 |
| Manuscript | `Manuscript/Revised Manuscript_2-1-26.docx` | Text with figure legends at end |

## Phase 1: Extraction

### 1.1 Figure Extraction (PowerPoint → PNG)

Uses `python-pptx` to extract embedded images from each slide.

```python
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE

def extract_figures(pptx_path, output_dir):
    prs = Presentation(pptx_path)
    for slide_idx, slide in enumerate(prs.slides):
        for shape in slide.shapes:
            if shape.shape_type == MSO_SHAPE_TYPE.PICTURE:
                # Save image blob
```

### 1.2 Legend Extraction (DOCX → JSON)

Parses figure legends from end of manuscript (paragraphs starting with "Fig." or "Extended Data Fig.").

Output format:
```json
{
  "Fig. 2": {
    "title": "Accelerated tumorigenesis in older rats...",
    "panels": {
      "A": "Oncoplot for three tumors from younger rats...",
      "B": "COSMIC signature analysis shows enrichment..."
    },
    "type": "COMPUTATIONAL",
    "analysis": ["rat_wes"]
  }
}
```

### 1.3 Claim Extraction

Parses quantitative claims from manuscript text:
- p-values: `p = 0.0033`, `p < 0.05`, `FDR < 0.05`
- Sample sizes: `n = 115`, `n = 89`
- Correlations: `Pearson r = 0.85`
- Percentages, fold changes

## Phase 2: Comparison

### Tier 1: Numerical Validation (Programmatic)

Compares extracted claims to code outputs:

| Status | Criteria |
|--------|----------|
| EXACT | Values match within tolerance (rtol=1e-3) |
| CLOSE | Same direction, within 10% |
| DIFFERS | Different values or direction |
| MISSING | Output not found |

### Tier 2: Visual Validation (Azure OpenAI GPT-5)

Uses Azure Foundry API for figure comparison:

```python
AZURE_OPENAI_ENDPOINT = 'https://alc37-meu8zr6o-eastus2.openai.azure.com/'
AZURE_OPENAI_DEPLOYMENT = 'gpt-5-chat'
AZURE_OPENAI_API_VERSION = '2024-12-01-preview'
AZURE_OPENAI_ENV_PATH = '/ihome/alee/alc376/.azure_openai_api_key'
```

Prompt structure:
1. Provide figure description from legend
2. Send both images (manuscript + regenerated)
3. Request structured JSON response with match level, differences, extracted values

Output format:
```json
{
  "figure_id": "Fig. 2A",
  "match_level": "EXACT|SIMILAR|DIFFERS",
  "confidence": 0.95,
  "data_patterns_match": true,
  "numerical_values_match": true,
  "differences": [],
  "extracted_values": {"SBS4": 0.45},
  "recommendation": "PASS|REVIEW|FAIL"
}
```

### Tier 3: Claim Validation

Maps manuscript claims to specific code outputs:
- Load relevant data files (CSV, RDS)
- Verify claimed statistics match computed values
- Check direction of effects (positive/negative correlation, up/down regulation)

## Phase 3: Report Generation

### Validation Report Structure

```markdown
# Manuscript Validation Report

## Executive Summary
| Tier | Total | PASS | REVIEW | FAIL | PENDING |

## Tier 1: Numerical Validation
[Per-claim comparison table]

## Tier 2: Visual Validation
[Per-figure comparison with GPT-5 analysis]

## Tier 3: Claim Validation
[Claim-to-data mapping results]

## Figures Requiring Review
[Flagged items needing manual verification]

## Appendix: Wet Lab Figures (Out of Scope)
[List of excluded figures]
```

### Status Definitions

| Status | Meaning | Action Required |
|--------|---------|-----------------|
| PASS | Exact or near-exact match | None |
| REVIEW | Minor differences, same conclusions | Manual verification |
| FAIL | Significant differences | Investigate discrepancy |
| PENDING | Output not yet generated | Run pipeline first |

## Scripts to Create

| Script | Purpose |
|--------|---------|
| `01_extract_figures.py` | PowerPoint → PNG extraction |
| `02_extract_claims.py` | DOCX → legends.json, claims.json |
| `02b_check_outputs.py` | Check which analysis outputs exist |
| `03_numerical_comparison.py` | CSV/RDS programmatic comparison |
| `04_visual_comparison.py` | Azure OpenAI GPT-5 vision comparison |
| `05_claim_verification.py` | Text claim → data verification |
| `06_generate_report.py` | Generate markdown report |
| `run_validation.sh` | Master orchestration script |

## Dependencies

Already available:
- `python-pptx` (PowerPoint parsing)
- `python-docx` (Word parsing)
- `openai` (Azure OpenAI client)
- `pandas`, `numpy` (data comparison)

## Notes

- Azure OpenAI GPT-5 vision tested and working
- Figure legends are at end of DOCX (paragraphs 400+)
- Main PowerPoint has 8 slides; Supplemental has 19 slides
- Some figures have multiple panels requiring separate extraction
