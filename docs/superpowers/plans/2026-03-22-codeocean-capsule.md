# ERpBRCA CodeOcean Capsule Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a CodeOcean capsule for the ERpBRCA_OlderWomen Nature Aging publication with full computational reproducibility from count matrices and a fast figure-regeneration mode.

**Architecture:** Section-isolated pipeline with 8 independent analysis modules, centralized config for path resolution with 3-tier fallback (intermediate → precomputed → figure_regen), dual-language R+Python Docker container, and a sync script to keep the capsule in sync with the main repo.

**Tech Stack:** R 4.3.3, Python 3.10, Bioconductor Docker, bash orchestration

**Spec:** `docs/superpowers/specs/2026-03-22-codeocean-capsule-design.md`

**Reference capsule:** `../../CITEgeistNeilAnalysis/CITEgeist_CodeOcean/`

---

## File Structure

### New files (in `../ERpBRCA_CodeOcean/`)

| File | Responsibility |
|------|---------------|
| `environment/Dockerfile` | bioconductor base + miniconda Python 3.10 |
| `environment/postInstall` | R + Python package installation |
| `code/config/paths.R` | Centralized path resolution with `resolve_input()` fallback |
| `code/config/paths.py` | Python equivalent of paths.R |
| `code/config/figure_config.R` | Shared R plot styling |
| `code/config/figure_config.py` | Shared Python plot styling |
| `code/run_all.sh` | Full computation orchestrator (trap-and-continue) |
| `code/run_codeocean.sh` | Figure-regeneration mode |
| `code/run_section.sh` | Single-section runner |
| `code/{01..07}_*/run.sh` | Per-section orchestrators (8 total) |
| `README.md` | Capsule documentation |
| `data/README.md` | Data provenance + checksums |
| `.gitignore` | Exclude data/, results/, __pycache__ |

### New files (in main repo `ERpBRCA_OlderWomen/`)

| File | Responsibility |
|------|---------------|
| `scripts/sync_to_capsule.sh` | One-directional code sync main→capsule |
| `scripts/prepare_data_asset.sh` | Package data asset (tier 1 or tier 2) |

### Files synced from main repo (not created, copied by sync script)

All analysis scripts per the inclusion lists in the spec. These are NOT modified — only config files are swapped.

---

## Task 1: Scaffold Capsule Directory

**Files:**
- Create: `../ERpBRCA_CodeOcean/` (entire directory tree)
- Create: `../ERpBRCA_CodeOcean/.gitignore`

- [ ] **Step 1: Create capsule directory structure**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview

mkdir -p ERpBRCA_CodeOcean/{environment,code/{config,01_human_bulk_rnaseq/mica,02_rat_snrnaseq,03_comparison,03_rat_wes,04_human_scrnaseq,05_rat_bulk_rnaseq,06_spatial_biopsies,07_organoid_single_cell},data,results/{figures,tables,intermediate,logs}}
```

- [ ] **Step 2: Create .gitignore**

```
data/
results/
__pycache__/
*.pyc
.ipynb_checkpoints/
*.Rhistory
```

- [ ] **Step 3: Initialize git repo**

```bash
cd ERpBRCA_CodeOcean && git init && git add .gitignore && git commit -m "chore: scaffold CodeOcean capsule directory"
```

---

## Task 2: Environment — Dockerfile and postInstall

**Files:**
- Create: `../ERpBRCA_CodeOcean/environment/Dockerfile`
- Create: `../ERpBRCA_CodeOcean/environment/postInstall`

- [ ] **Step 1: Write Dockerfile**

Base image: `bioconductor/bioconductor_docker:3.19` (R 4.3.3, Ubuntu 22.04). Add miniconda for Python 3.10.

```dockerfile
FROM bioconductor/bioconductor_docker:3.19

# Install miniconda for Python 3.10
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH
RUN conda install -y python=3.10 && conda clean -afy
```

- [ ] **Step 2: Write postInstall**

This runs after Docker build. Installs all R and Python packages.

```bash
#!/usr/bin/env bash
set -e

echo "=== Installing R packages ==="

# Bioconductor packages
Rscript -e '
BiocManager::install(c(
  "DESeq2", "GSVA", "PROGENy", "SingleR", "SingleCellExperiment",
  "msigdbr", "qusage", "speckle"
), ask=FALSE, update=FALSE)
'

# CRAN packages
Rscript -e '
install.packages(c(
  "Seurat", "harmony", "DoubletFinder", "genefu",
  "patchwork", "pheatmap", "ggplot2", "data.table",
  "lme4", "dittoseq", "tidyverse"
), repos="https://cloud.r-project.org", Ncpus=4)
'

# CellChat (requires devtools)
Rscript -e '
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
devtools::install_github("jinworks/CellChat")
'

echo "=== Installing Python packages ==="
pip install --no-cache-dir \
  scanpy anndata squidpy gseapy pydeseq2 decoupler \
  scrublet statsmodels adjustText psutil scikit-misc \
  matplotlib seaborn

echo "=== Installation complete ==="
```

**Note:** CellPhoneDB, commot, and MICA are deferred. During Task 8 (dependency audit), we check whether Sec 06 scripts import commot directly. If so, precomputed results are used; the scripts get a guard clause.

- [ ] **Step 3: Commit**

```bash
git add environment/ && git commit -m "feat: add Dockerfile and postInstall for R 4.3.3 + Python 3.10"
```

---

## Task 3: Centralized Config — paths.R and paths.py

**Files:**
- Create: `../ERpBRCA_CodeOcean/code/config/paths.R`
- Create: `../ERpBRCA_CodeOcean/code/config/paths.py`

- [ ] **Step 1: Write paths.R**

This is the core path resolution layer. Every R script sources this at the top.

```r
# code/config/paths.R — Centralized path configuration for CodeOcean capsule
#
# Usage: source(file.path(dirname(sys.frame(1)$ofile), "..", "config", "paths.R"))
# Or:    source("code/config/paths.R")  # from run.sh working directory

# Detect capsule root (parent of code/)
.find_root <- function() {
  # Try script directory first
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("--file=", "", file_arg)))
    # Walk up to find code/ parent
    candidate <- script_dir
    for (i in 1:5) {
      if (file.exists(file.path(candidate, "code", "config", "paths.R"))) {
        return(candidate)
      }
      candidate <- dirname(candidate)
    }
  }
  # Fallback: working directory
  return(getwd())
}

ROOT_DIR   <- .find_root()
CODE_DIR   <- file.path(ROOT_DIR, "code")
DATA_DIR   <- file.path(ROOT_DIR, "data")
RESULTS_DIR <- file.path(ROOT_DIR, "results")

# CodeOcean mount overrides: /code, /data, /results are separate mounts
if (dir.exists("/code") && !dir.exists(CODE_DIR)) CODE_DIR <- "/code"
if (dir.exists("/data") && !dir.exists(DATA_DIR)) DATA_DIR <- "/data"
if (dir.exists("/results") && !dir.exists(RESULTS_DIR)) RESULTS_DIR <- "/results"

# Resolve an input file with 3-tier fallback:
#   1. results/intermediate/  (from current run_all.sh)
#   2. data/precomputed/      (Tier 2 large objects)
#   3. data/figure_regen/     (Tier 1 lightweight CSVs)
resolve_input <- function(section, filename) {
  candidates <- c(
    file.path(RESULTS_DIR, "intermediate", paste0(section, "_outputs"), filename),
    file.path(DATA_DIR, "precomputed", paste0(section, "_outputs"), filename),
    file.path(DATA_DIR, "figure_regen", paste0(section, "_outputs"), filename)
  )
  for (path in candidates) {
    if (file.exists(path)) {
      cat("  Input resolved:", path, "\n")
      return(path)
    }
  }
  stop(paste("No input found for", section, "/", filename,
             "\nSearched:", paste(candidates, collapse="\n  ")))
}

# Save an intermediate output (always to results/intermediate/)
save_intermediate <- function(section, filename) {
  dir <- file.path(RESULTS_DIR, "intermediate", paste0(section, "_outputs"))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  file.path(dir, filename)
}

# Save a figure (to results/figures/<section>/)
save_figure <- function(section, filename) {
  dir <- file.path(RESULTS_DIR, "figures", section)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  file.path(dir, filename)
}

# Section 03_comparison special case: reads from data/03_comparison/
COMPARISON_ORIGINAL_DIR  <- file.path(DATA_DIR, "03_comparison", "original")
COMPARISON_CORRECTED_DIR <- file.path(DATA_DIR, "03_comparison", "corrected")

cat("=== Capsule paths initialized ===\n")
cat("  ROOT:", ROOT_DIR, "\n")
cat("  DATA:", DATA_DIR, "\n")
cat("  RESULTS:", RESULTS_DIR, "\n\n")
```

- [ ] **Step 2: Write paths.py**

```python
"""code/config/paths.py — Centralized path configuration for CodeOcean capsule.

Usage:
    import sys; sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "config"))
    from paths import ROOT_DIR, DATA_DIR, RESULTS_DIR, resolve_input, save_intermediate, save_figure
"""
from pathlib import Path

# Detect capsule root
_config_dir = Path(__file__).resolve().parent  # code/config/
CODE_DIR = _config_dir.parent                   # code/
ROOT_DIR = CODE_DIR.parent                      # ERpBRCA_CodeOcean/

DATA_DIR = ROOT_DIR / "data"
RESULTS_DIR = ROOT_DIR / "results"

# CodeOcean override: /data and /results may be mounted
if Path("/data").exists() and not DATA_DIR.exists():
    DATA_DIR = Path("/data")
if Path("/results").exists() and not RESULTS_DIR.exists():
    RESULTS_DIR = Path("/results")


def resolve_input(section: str, filename: str) -> Path:
    """Resolve input with 3-tier fallback."""
    candidates = [
        RESULTS_DIR / "intermediate" / f"{section}_outputs" / filename,
        DATA_DIR / "precomputed" / f"{section}_outputs" / filename,
        DATA_DIR / "figure_regen" / f"{section}_outputs" / filename,
    ]
    for path in candidates:
        if path.exists():
            print(f"  Input resolved: {path}")
            return path
    searched = "\n  ".join(str(c) for c in candidates)
    raise FileNotFoundError(f"No input found for {section}/{filename}\nSearched:\n  {searched}")


def save_intermediate(section: str, filename: str) -> Path:
    """Return path for saving an intermediate output."""
    out_dir = RESULTS_DIR / "intermediate" / f"{section}_outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / filename


def save_figure(section: str, filename: str) -> Path:
    """Return path for saving a figure."""
    fig_dir = RESULTS_DIR / "figures" / section
    fig_dir.mkdir(parents=True, exist_ok=True)
    return fig_dir / filename


# Section 03_comparison special case
COMPARISON_ORIGINAL_DIR = DATA_DIR / "03_comparison" / "original"
COMPARISON_CORRECTED_DIR = DATA_DIR / "03_comparison" / "corrected"

print(f"=== Capsule paths initialized ===")
print(f"  ROOT: {ROOT_DIR}")
print(f"  DATA: {DATA_DIR}")
print(f"  RESULTS: {RESULTS_DIR}")
```

- [ ] **Step 3: Commit**

```bash
git add code/config/paths.R code/config/paths.py
git commit -m "feat: add centralized path config with 3-tier input resolution"
```

---

## Task 4: Figure Config Files

**Files:**
- Create: `../ERpBRCA_CodeOcean/code/config/figure_config.R`
- Create: `../ERpBRCA_CodeOcean/code/config/figure_config.py`

- [ ] **Step 1: Copy figure configs from main repo**

Check existing figure config patterns in the main repo:

```bash
# Find existing figure config files
cat ERpBRCA_OlderWomen/analysis/06_spatial_biopsies/figure_config.py
cat ERpBRCA_OlderWomen/analysis/07_organoid_single_cell/_figure_config.py
```

Adapt these into shared capsule versions. For R sections, check if any common theme/palette is used across scripts 01-05.

- [ ] **Step 2: Write figure_config.R**

Extract common ggplot2 theme and color palettes from existing R scripts. Minimal — only shared styling.

- [ ] **Step 3: Write figure_config.py**

Merge `figure_config.py` and `_figure_config.py` patterns from sections 06/07 into one shared module.

- [ ] **Step 4: Commit**

```bash
git add code/config/figure_config.* && git commit -m "feat: add shared figure config for R and Python"
```

---

## Task 5: Orchestrator Scripts — run_all.sh, run_codeocean.sh, run_section.sh

**Files:**
- Create: `../ERpBRCA_CodeOcean/code/run_all.sh`
- Create: `../ERpBRCA_CodeOcean/code/run_codeocean.sh`
- Create: `../ERpBRCA_CodeOcean/code/run_section.sh`

- [ ] **Step 1: Write run_all.sh**

Trap-and-continue pattern. Does NOT use `set -e`. Logs per-section status.

```bash
#!/bin/bash
# run_all.sh — Full computation from count matrices (~4-6 hours)
# Runs all 8 sections sequentially. Continues on section failure.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$ROOT_DIR/results/logs"
mkdir -p "$LOG_DIR"

SECTIONS=(
  01_human_bulk_rnaseq
  02_rat_snrnaseq
  03_comparison
  03_rat_wes
  04_human_scrnaseq
  05_rat_bulk_rnaseq
  06_spatial_biopsies
  07_organoid_single_cell
)

FAILED=0
RESULTS=()

echo "=============================================="
echo "ERpBRCA Analysis Pipeline — Full Computation"
echo "=============================================="
echo "Start time: $(date)"
echo ""

for SECTION in "${SECTIONS[@]}"; do
  echo "--- [$SECTION] Starting at $(date) ---"
  RUN_SH="$SCRIPT_DIR/$SECTION/run.sh"

  if [ ! -f "$RUN_SH" ]; then
    echo "  SKIP: $RUN_SH not found"
    RESULTS+=("$SECTION: SKIPPED")
    continue
  fi

  if bash "$RUN_SH" > >(tee "$LOG_DIR/${SECTION}.log") 2>&1; then
    echo "  [$SECTION] PASSED"
    RESULTS+=("$SECTION: PASSED")
  else
    echo "  [$SECTION] FAILED (see $LOG_DIR/${SECTION}.log)"
    RESULTS+=("$SECTION: FAILED")
    FAILED=$((FAILED + 1))
  fi
  echo ""
done

echo "=============================================="
echo "Pipeline Summary — $(date)"
echo "=============================================="
for R in "${RESULTS[@]}"; do echo "  $R"; done
echo ""
echo "Results in: $ROOT_DIR/results/"
du -sh "$ROOT_DIR/results/" 2>/dev/null || true

exit $FAILED
```

- [ ] **Step 2: Write run_codeocean.sh**

Figure-regeneration mode. Sets `CODEOCEAN_FIGURE_REGEN=1` env var that section `run.sh` scripts check.

```bash
#!/bin/bash
# run_codeocean.sh — Figure regeneration from precomputed data (~15-30 min)
# Reads from data/precomputed/ or data/figure_regen/ and regenerates all figures.

export CODEOCEAN_FIGURE_REGEN=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "=============================================="
echo "ERpBRCA Analysis — Figure Regeneration Mode"
echo "=============================================="
echo ""

exec bash "$SCRIPT_DIR/run_all.sh"
```

- [ ] **Step 3: Write run_section.sh**

```bash
#!/bin/bash
# run_section.sh — Run a single analysis section
# Usage: ./run_section.sh 02_rat_snrnaseq
#        ./run_section.sh 03_rat_wes

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -z "$1" ]; then
  echo "Usage: $0 <section_name>"
  echo ""
  echo "Available sections:"
  for d in "$SCRIPT_DIR"/*/; do
    [ -f "$d/run.sh" ] && echo "  $(basename "$d")"
  done
  exit 1
fi

SECTION="$1"
RUN_SH="$SCRIPT_DIR/$SECTION/run.sh"

if [ ! -f "$RUN_SH" ]; then
  echo "Error: $RUN_SH not found"
  exit 1
fi

echo "Running section: $SECTION"
exec bash "$RUN_SH"
```

- [ ] **Step 4: Make scripts executable and commit**

```bash
chmod +x code/run_all.sh code/run_codeocean.sh code/run_section.sh
git add code/run_all.sh code/run_codeocean.sh code/run_section.sh
git commit -m "feat: add orchestrator scripts (run_all, run_codeocean, run_section)"
```

---

## Task 6: Per-Section run.sh Files

**Files:**
- Create: `../ERpBRCA_CodeOcean/code/{section}/run.sh` (8 files)

Each `run.sh` has a header comment documenting language, runtime, inputs, outputs, and manuscript figures. Uses `set -e`. Checks `CODEOCEAN_FIGURE_REGEN` to decide whether to run full computation or figure-only scripts.

- [ ] **Step 1: Write run.sh for 01_human_bulk_rnaseq**

```bash
#!/bin/bash
# Section 01: Human Bulk RNA-seq Analysis
# Language: R
# Runtime: ~20 min (full), ~5 min (figure-regen)
# Inputs: data/human_bulk_rnaseq/ (counts, TPM, metadata)
# Outputs: results/figures/01_human_bulk_rnaseq/, results/intermediate/01_outputs/
# Manuscript: Figures 1, 4, 5

set -e
SECTION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[01_human_bulk_rnaseq] Starting..."

if [ "${CODEOCEAN_FIGURE_REGEN:-0}" = "1" ]; then
  echo "  Figure regeneration mode — running visualization only"
  Rscript "$SECTION_DIR/05_visualize.R"
  Rscript "$SECTION_DIR/mica/04_fig4_estrogen.R"
  Rscript "$SECTION_DIR/mica/05_fig5_immune.R"
else
  Rscript "$SECTION_DIR/01_preprocess.R"
  Rscript "$SECTION_DIR/02_run_gsva.R"
  Rscript "$SECTION_DIR/03_run_progeny.R"
  Rscript "$SECTION_DIR/04_correlations.R"
  Rscript "$SECTION_DIR/05_visualize.R"
  echo "  Running MICA subpipeline..."
  Rscript "$SECTION_DIR/mica/01_prep_data.R"
  Rscript "$SECTION_DIR/mica/02_run_gsva.R"
  Rscript "$SECTION_DIR/mica/03_run_mica.R"
  Rscript "$SECTION_DIR/mica/04_fig4_estrogen.R"
  Rscript "$SECTION_DIR/mica/05_fig5_immune.R"
  Rscript "$SECTION_DIR/mica/07_statistical_audit.R"
fi

echo "[01_human_bulk_rnaseq] Complete."
```

- [ ] **Step 2: Write run.sh for 02_rat_snrnaseq**

Same pattern. Full mode runs all 9 scripts. Figure-regen mode: check each script to determine which produce figures vs compute. The implementer must `grep` for `ggsave\|pdf(\|png(\|DimPlot\|FeaturePlot\|VlnPlot\|DoHeatmap` in each script. Scripts that ONLY produce figures (no `saveRDS`) go in figure-regen. Scripts that compute + plot go in full-only (their plots are regenerated from precomputed intermediates). If a script does both (compute + save figure), split the figure-saving into the figure-regen branch using `resolve_input()` to load the precomputed object.

- [ ] **Step 3: Write run.sh for 03_comparison**

```bash
#!/bin/bash
# Section 03_comparison: Cross-Species Comparison
# Language: R
# Runtime: ~2 min
# Inputs: data/03_comparison/{original,corrected}/
# Outputs: results/figures/03_comparison/

set -e
SECTION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[03_comparison] Starting..."
Rscript "$SECTION_DIR/01_generate_comparison_report.R"
echo "[03_comparison] Complete."
```

No figure-regen distinction needed — this script is fast and always reads from data/.

- [ ] **Step 4: Write run.sh for 03_rat_wes**

```bash
#!/bin/bash
# Section 03_rat_wes: Rat Whole Exome Sequencing
# Language: Python
# Runtime: ~10 min
# Inputs: data/rat_wes/ (VEP-annotated outputs)
# Outputs: results/figures/03_rat_wes/

set -e
SECTION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[03_rat_wes] Starting..."

if [ "${CODEOCEAN_FIGURE_REGEN:-0}" = "1" ]; then
  python "$SECTION_DIR/03_generate_oncoplot.py"
else
  python "$SECTION_DIR/01_parse_vep.py"
  python "$SECTION_DIR/02_cosmic_signatures.py"
  python "$SECTION_DIR/03_generate_oncoplot.py"
fi

echo "[03_rat_wes] Complete."
```

- [ ] **Step 5: Write run.sh for 04_human_scrnaseq**

Largest section (19 scripts). Figure-regen mode runs only visualization scripts (04, 05, 06, 07, 12, 15, 16, 17). Full mode runs all 19.

```bash
#!/bin/bash
# Section 04: Human scRNA-seq Analysis (Xu et al. 2024)
# Language: R + Python
# Runtime: ~2-4 hours (full), ~15 min (figure-regen)
# Inputs: data/human_scrnaseq/ (Xu et al. count matrix + metadata)
# Outputs: results/figures/04_human_scrnaseq/

set -e
SECTION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[04_human_scrnaseq] Starting..."

if [ "${CODEOCEAN_FIGURE_REGEN:-0}" = "1" ]; then
  echo "  Figure regeneration mode — plotting from precomputed intermediates"
  # Only scripts that produce publication figures from precomputed data.
  # Scripts 00b, 01, 03, 06-14, 16 are computation-only or produce intermediates.
  # The implementer MUST audit each script: grep for ggsave/pdf/png calls
  # and verify it can load from resolve_input() before adding here.
  Rscript "$SECTION_DIR/04_cell_fractions.R"
  Rscript "$SECTION_DIR/05_gene_expression_violin.R"
  Rscript "$SECTION_DIR/15_multicelltype_pathway.R"
  Rscript "$SECTION_DIR/17_cellphonedb_dotplot.R"
  # NOTE: This list is a starting point. During Task 11 (path adaptation),
  # verify which scripts produce manuscript figures and add them here.
else
  Rscript "$SECTION_DIR/00b_preprocess_seurat.R"
  Rscript "$SECTION_DIR/01_load_subset_data.R"
  Rscript "$SECTION_DIR/03_cell_type_annotation.R"
  Rscript "$SECTION_DIR/04_cell_fractions.R"
  Rscript "$SECTION_DIR/05_gene_expression_violin.R"
  Rscript "$SECTION_DIR/06_run_gsva.R"
  Rscript "$SECTION_DIR/07_run_progeny.R"
  Rscript "$SECTION_DIR/08_run_wcsea.R"
  Rscript "$SECTION_DIR/09_cellphonedb_prep.R"
  Rscript "$SECTION_DIR/10_load_macrophage_seurat.R"
  Rscript "$SECTION_DIR/11_macrophage_deg_analysis.R"
  Rscript "$SECTION_DIR/12_macrophage_pathway_enrichment.R"
  python  "$SECTION_DIR/12_enrichr_pathway_analysis.py"
  Rscript "$SECTION_DIR/13_macrophage_cellphonedb.R"
  Rscript "$SECTION_DIR/14_multicelltype_deg.R"
  Rscript "$SECTION_DIR/15_multicelltype_pathway.R"
  Rscript "$SECTION_DIR/16_cellchat_analysis.R"
  python  "$SECTION_DIR/16_run_cellphonedb.py"
  Rscript "$SECTION_DIR/17_cellphonedb_dotplot.R"
fi

echo "[04_human_scrnaseq] Complete."
```

- [ ] **Step 6: Write run.sh for 05_rat_bulk_rnaseq**

```bash
#!/bin/bash
# Section 05: Rat Bulk RNA-seq (post-alignment)
# Language: R
# Runtime: ~15 min
# Inputs: data/rat_bulk_rnaseq/ (htseq count matrices)

set -e
SECTION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[05_rat_bulk_rnaseq] Starting..."

Rscript "$SECTION_DIR/03_deseq2.R"
Rscript "$SECTION_DIR/04_pam50_subtyping.R"
Rscript "$SECTION_DIR/05_validate_vs_rahul.R"

echo "[05_rat_bulk_rnaseq] Complete."
```

In figure-regen mode, `03_deseq2.R` must load precomputed DESeq2 results via `resolve_input()` rather than re-fitting. Add a figure-regen guard:

```bash
if [ "${CODEOCEAN_FIGURE_REGEN:-0}" = "1" ]; then
  echo "  Figure regeneration mode — plotting from precomputed results"
  Rscript "$SECTION_DIR/04_pam50_subtyping.R"
  Rscript "$SECTION_DIR/05_validate_vs_rahul.R"
else
  Rscript "$SECTION_DIR/03_deseq2.R"
  Rscript "$SECTION_DIR/04_pam50_subtyping.R"
  Rscript "$SECTION_DIR/05_validate_vs_rahul.R"
fi
```

Scripts 04 and 05 must use `resolve_input()` to load DESeq2 results either from `results/intermediate/` or `data/figure_regen/`.

- [ ] **Step 7: Write run.sh for 06_spatial_biopsies**

```bash
#!/bin/bash
# Section 06: Spatial Biopsies (HCC22-088)
# Language: Python
# Runtime: ~10 min
# Inputs: data/spatial_biopsies/

set -e
SECTION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[06_spatial_biopsies] Starting..."

# Both scripts read from data/ and produce figures — no heavy computation
# In figure-regen mode, they read from data/figure_regen/ if spatial data unavailable
python "$SECTION_DIR/01_immune_secretion.py"
python "$SECTION_DIR/02_immune_pathways.py"

echo "[06_spatial_biopsies] Complete."
```

- [ ] **Step 8: Write run.sh for 07_organoid_single_cell**

```bash
#!/bin/bash
# Section 07: Organoid Single Cell (PDO-296, post-Cell Ranger)
# Language: Python
# Runtime: ~30 min (full), ~5 min (figure-regen)
# Inputs: data/organoid_single_cell/ (filtered matrices)
# Note: Uses seed=42 (not 12345) for consistency with original results

set -e
SECTION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[07_organoid_single_cell] Starting..."

if [ "${CODEOCEAN_FIGURE_REGEN:-0}" = "1" ]; then
  echo "  Figure regeneration mode"
  python "$SECTION_DIR/08_inhibitor_mechanism_exploration.py"
  python "$SECTION_DIR/09_pathway_enrichment_analysis.py"
  python "$SECTION_DIR/10_heterogeneity_analysis.py"
  python "$SECTION_DIR/11_proliferation_analysis.py"
else
  python "$SECTION_DIR/02_qc.py"
  python "$SECTION_DIR/03_preprocess.py"
  python "$SECTION_DIR/04_pseudobulk_de.py"
  python "$SECTION_DIR/05_pathways.py"
  python "$SECTION_DIR/06_cell_cycle.py"
  python "$SECTION_DIR/07_single_cell_pathways.py"
  python "$SECTION_DIR/08_inhibitor_mechanism_exploration.py"
  python "$SECTION_DIR/09_pathway_enrichment_analysis.py"
  python "$SECTION_DIR/10_heterogeneity_analysis.py"
  python "$SECTION_DIR/11_proliferation_analysis.py"
fi

echo "[07_organoid_single_cell] Complete."
```

- [ ] **Step 9: Make all run.sh executable and commit**

```bash
find code/ -name "run.sh" -exec chmod +x {} \;
git add code/*/run.sh
git commit -m "feat: add per-section run.sh orchestrators with figure-regen mode"
```

---

## Task 7: Sync Script

**Files:**
- Create: `ERpBRCA_OlderWomen/scripts/sync_to_capsule.sh`

- [ ] **Step 1: Write sync_to_capsule.sh**

Key behaviors:
- Reads inclusion/exclusion lists (hardcoded arrays matching the spec)
- Copies scripts from `analysis/XX_*/` → `../ERpBRCA_CodeOcean/code/XX_*/`
- Skips SLURM wrappers, alignment scripts, `__pycache__/`, `.pipeline_markers/`
- Does NOT copy `run.sh` (those are capsule-specific, not from main repo)
- Does NOT copy config files (capsule has its own `paths.R`/`paths.py`)
- Generates `MANIFEST.md` with source commit hash per file
- `--dry-run`, `--section`, `--diff` flags
- `--target` to override capsule directory

The script defines per-section inclusion arrays:

```bash
declare -A INCLUDE
INCLUDE[01_human_bulk_rnaseq]="01_preprocess.R 02_run_gsva.R 03_run_progeny.R 04_correlations.R 05_visualize.R"
INCLUDE[01_human_bulk_rnaseq/mica]="01_prep_data.R 02_run_gsva.R 03_run_mica.R 04_fig4_estrogen.R 05_fig5_immune.R 07_statistical_audit.R"
INCLUDE[02_rat_snrnaseq]="01_qc_filter.R 02_normalize_integrate.R 03_cluster_annotate.R 03b_sctype_annotate.R 03c_sctype_cluster.R 08_differential_expression.R 09_differential_abundance.R method_comparison.R spot_check_markers.R"
INCLUDE[03_comparison]="01_generate_comparison_report.R"
INCLUDE[03_rat_wes]="01_parse_vep.py 02_cosmic_signatures.py 03_generate_oncoplot.py"
INCLUDE[04_human_scrnaseq]="00b_preprocess_seurat.R 01_load_subset_data.R 03_cell_type_annotation.R 04_cell_fractions.R 05_gene_expression_violin.R 06_run_gsva.R 07_run_progeny.R 08_run_wcsea.R 09_cellphonedb_prep.R 10_load_macrophage_seurat.R 11_macrophage_deg_analysis.R 12_macrophage_pathway_enrichment.R 12_enrichr_pathway_analysis.py 13_macrophage_cellphonedb.R 14_multicelltype_deg.R 15_multicelltype_pathway.R 16_cellchat_analysis.R 16_run_cellphonedb.py 17_cellphonedb_dotplot.R"
INCLUDE[05_rat_bulk_rnaseq]="03_deseq2.R 04_pam50_subtyping.R 05_validate_vs_rahul.R"
INCLUDE[06_spatial_biopsies]="01_immune_secretion.py 02_immune_pathways.py figure_config.py utils.py"
INCLUDE[07_organoid_single_cell]="02_qc.py 03_preprocess.py 04_pseudobulk_de.py 05_pathways.py 06_cell_cycle.py 07_single_cell_pathways.py 08_inhibitor_mechanism_exploration.py 09_pathway_enrichment_analysis.py 10_heterogeneity_analysis.py 11_proliferation_analysis.py _config.py _figure_config.py"
```

After copying, for each section, inject a `source("../../config/paths.R")` line at the top of R scripts that currently use hardcoded HPC paths (checked via grep).

- [ ] **Step 2: Test dry-run**

```bash
./scripts/sync_to_capsule.sh --dry-run
```

Verify: correct files listed, no SLURM wrappers, no alignment scripts.

- [ ] **Step 3: Run actual sync**

```bash
./scripts/sync_to_capsule.sh
```

Verify: files appear in `../ERpBRCA_CodeOcean/code/`, MANIFEST.md generated.

- [ ] **Step 4: Commit sync script in main repo**

```bash
cd ERpBRCA_OlderWomen
git add scripts/sync_to_capsule.sh
git commit -m "feat: add sync_to_capsule.sh for one-directional code sync"
```

---

## Task 8: Dependency Audit — Resolve Decision Points

**Files:** May modify `environment/postInstall`, section `run.sh` files

This task resolves the 3 open dependency questions from the spec.

- [ ] **Step 1: Check if Sec 06 scripts import commot**

```bash
grep -n "import commot\|from commot" ERpBRCA_OlderWomen/analysis/06_spatial_biopsies/*.py
```

If they do: add a guard clause that checks for commot availability and falls back to precomputed.
If they don't: no action needed.

- [ ] **Step 2: Check MICA availability**

```bash
Rscript -e 'install.packages("MICA"); library(MICA)'
# Or check BiocManager
Rscript -e 'BiocManager::available("MICA")'
```

If available: add to postInstall.
If not: ensure `mica/03_run_mica.R` has a precomputed fallback (reads MICA output from data asset), and figure scripts (04, 05) can run without re-running MICA.

- [ ] **Step 3: Check CellPhoneDB v5 install complexity**

```bash
pip install cellphonedb 2>&1 | tail -20
```

If clean: add to postInstall.
If not: ensure scripts 09, 13, 16 have precomputed fallbacks, and script 17 (visualization) can read saved results.

- [ ] **Step 4: Update postInstall and run.sh files based on findings**

- [ ] **Step 5: Commit changes**

```bash
cd ../ERpBRCA_CodeOcean
git add -A && git commit -m "feat: resolve dependency decision points (MICA, CellPhoneDB, commot)"
```

---

## Task 9: Data Packaging Script

**Files:**
- Create: `ERpBRCA_OlderWomen/scripts/prepare_data_asset.sh`

- [ ] **Step 1: Write prepare_data_asset.sh**

Takes `--tier 1` or `--tier 2`. Key logic:

```bash
# Tier 1 data sources (from ERpBRCA_OlderWomen/):
#   data/human_bulk_rnaseq/raw/         → human_bulk_rnaseq/
#   data/rat_snrnaseq/raw/              → rat_snrnaseq/
#   data/rat_wes/raw/                   → rat_wes/
#   /ix1/.../Xu_etal_.../matrix.mtx     → human_scrnaseq/ (gzip compressed)
#   data/rat_bulk_rnaseq/               → rat_bulk_rnaseq/
#   <spatial biopsy data>               → spatial_biopsies/
#   <post-CR matrices for Sec 07>       → organoid_single_cell/
#   data/gmt/                           → gmt/
#   results/original + results/corrected → 03_comparison/
#   <figure-regen CSVs>                 → figure_regen/
```

The script:
1. Creates a staging directory
2. Copies/compresses inputs per tier
3. For Tier 1: exports figure-regen CSVs from existing analysis outputs. This requires identifying which CSV/TSV files each section's plotting scripts actually read. Implementation approach: grep plotting scripts for `read.csv`/`read_csv`/`pd.read_csv`/`readRDS` calls and trace their file paths.
4. Generates `data/README.md` with SHA256 checksums
5. Reports total size
6. Validates all sections have required inputs

- [ ] **Step 2: Test tier 1 packaging**

```bash
./scripts/prepare_data_asset.sh --tier 1 --output /tmp/erp_data_tier1
du -sh /tmp/erp_data_tier1
```

Verify total size < 2.5 GB.

- [ ] **Step 3: Test tier 2 packaging** (if applicable)

```bash
./scripts/prepare_data_asset.sh --tier 2 --output /tmp/erp_data_tier2
du -sh /tmp/erp_data_tier2
```

- [ ] **Step 4: Commit**

```bash
cd ERpBRCA_OlderWomen
git add scripts/prepare_data_asset.sh
git commit -m "feat: add prepare_data_asset.sh for tier 1/2 CodeOcean data packaging"
```

---

## Task 10: README Files

**Files:**
- Create: `../ERpBRCA_CodeOcean/README.md`
- Create: `../ERpBRCA_CodeOcean/data/README.md` (template — checksums populated by prepare_data_asset.sh)

- [ ] **Step 1: Write capsule README.md**

Follow the structure from the spec:
- Title, authors, journal
- Quick Start (3 commands)
- Analysis Sections table (section, description, language, runtime, manuscript figures)
- Data Provenance table (dataset, GEO accession, citation, preprocessing)
- Excluded Steps table with exact commands/parameters
- Environment (R/Python versions, packages)
- Output Map (results path → manuscript figure)
- Seeding & Reproducibility

Reference the CITEgeist capsule README for formatting patterns but adapt for 8 sections and mixed R/Python.

- [ ] **Step 2: Write data/README.md template**

Provenance for each data subdirectory. Checksums section left as `<!-- Generated by prepare_data_asset.sh -->` placeholder.

- [ ] **Step 3: Commit**

```bash
cd ../ERpBRCA_CodeOcean
git add README.md data/README.md
git commit -m "docs: add capsule README and data provenance template"
```

---

## Task 11: Script Path Adaptation

**Files:** Modify synced analysis scripts in `../ERpBRCA_CodeOcean/code/`

After the sync (Task 7), the analysis scripts still contain HPC-specific paths. This task adapts them.

**Special cases to handle explicitly:**
- **Section 03_comparison:** Does NOT use `resolve_input()`. Must replace hardcoded `results/original/` and `results/corrected/` paths with `COMPARISON_ORIGINAL_DIR` and `COMPARISON_CORRECTED_DIR` from `paths.R`.
- **Section 07 `_config.py`:** The sync script DOES copy this file (it's in the inclusion list). It must be adapted to use capsule paths — replace HPC-specific data paths and the `standardize_treatments()` function's path references. The capsule's `paths.py` provides `DATA_DIR` etc.
- **Section 04 `00a_download_xu_geo.sh`:** Create a NEW capsule-specific download script at `code/04_human_scrnaseq/00a_download_xu_geo.sh`. This downloads the Xu et al. 2024 data from GEO if not present in the data asset. The GEO accession must be found from the paper's key resources table (doi:10.1016/j.xcrm.2024.101500).

- [ ] **Step 1: Identify HPC paths in synced scripts**

```bash
cd ../ERpBRCA_CodeOcean
grep -rn "/ix1/alee\|/ihome/alee\|/bgfs/" code/ --include="*.R" --include="*.py" | head -40
```

- [ ] **Step 2: For each R script, add paths.R source line**

Each R script needs to source the centralized config. Add at the top (after the shebang/header comments):

```r
# Source capsule path configuration
source(file.path(dirname(sys.frame(1)$ofile), "..", "config", "paths.R"))
```

Then replace hardcoded paths with `resolve_input()`, `save_intermediate()`, `save_figure()` calls.

**Strategy:** This is the most labor-intensive task. Approach per-section:
1. Read each script to find all `file.path()`, `readRDS()`, `read.csv()`, `ggsave()` calls
2. Replace input paths with `resolve_input("section", "filename")`
3. Replace output paths with `save_intermediate()` / `save_figure()`
4. Test that the script loads without error (syntax check only, not execution)

- [ ] **Step 3: For each Python script, add paths.py import**

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "config"))
from paths import resolve_input, save_intermediate, save_figure, DATA_DIR
```

- [ ] **Step 4: Create 00a_download_xu_geo.sh for Section 04**

Write a capsule-specific download script that fetches Xu et al. 2024 data from GEO if not already in `data/human_scrnaseq/`. Uses `wget` with the GEO FTP URL. The script is a fallback for Tier 1 if Xu data was not included in the data asset.

```bash
#!/bin/bash
# Download Xu et al. 2024 Primary Breast Tumor Atlas from GEO
# Only needed if data is not bundled in the data asset
set -e
DATA_DIR="${1:-../data/human_scrnaseq}"
if [ -f "$DATA_DIR/matrix.mtx" ] || [ -f "$DATA_DIR/matrix.mtx.gz" ]; then
  echo "Xu et al. data already present, skipping download"
  exit 0
fi
echo "Downloading Xu et al. 2024 data from GEO..."
mkdir -p "$DATA_DIR"
# TODO: Replace with actual GEO accession URL from paper key resources table
# wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/..." -O "$DATA_DIR/matrix.mtx.gz"
echo "ERROR: GEO accession URL not yet configured. See spec open question #5."
exit 1
```

- [ ] **Step 5: Adapt Section 03_comparison paths**

Replace `results/original/` and `results/corrected/` with `COMPARISON_ORIGINAL_DIR` and `COMPARISON_CORRECTED_DIR` in the script.

- [ ] **Step 6: Adapt Section 07 _config.py**

Replace HPC-specific paths in `_config.py` with imports from capsule `paths.py`. Keep `standardize_treatments()` function intact.

- [ ] **Step 7: Verify no HPC paths remain**

```bash
grep -rn "/ix1/\|/ihome/\|/bgfs/" code/ --include="*.R" --include="*.py"
```

Should return 0 results.

- [ ] **Step 8: Commit per section** (8 commits, one per section)

```bash
git add code/01_human_bulk_rnaseq/ && git commit -m "refactor(01): adapt paths for CodeOcean capsule"
# ... repeat for each section
```

---

## Task 12: End-to-End Smoke Test

**Files:** None (testing only)

- [ ] **Step 1: Verify capsule structure**

```bash
cd ../ERpBRCA_CodeOcean
find . -type f | grep -v ".git/" | sort
```

Confirm all expected files exist per spec directory tree.

- [ ] **Step 2: Verify no HPC paths**

```bash
grep -rn "/ix1/\|/ihome/\|/bgfs/" code/ --include="*.R" --include="*.py"
```

- [ ] **Step 3: Verify git repo size < 1 GB**

```bash
du -sh .git/
find . -size +100M -not -path "./.git/*"
```

No individual file > 100 MB.

- [ ] **Step 4: Syntax-check R scripts**

```bash
for f in $(find code/ -name "*.R"); do
  Rscript -e "parse('$f')" 2>&1 | grep -i error && echo "FAIL: $f"
done
```

- [ ] **Step 5: Syntax-check Python scripts**

```bash
for f in $(find code/ -name "*.py"); do
  python -m py_compile "$f" 2>&1 && echo "OK: $f" || echo "FAIL: $f"
done
```

- [ ] **Step 6: Test run_section.sh with no data (should fail gracefully)**

```bash
./code/run_section.sh 03_rat_wes 2>&1 | head -20
```

Should fail with "No input found" from `resolve_input()`, not a cryptic error.

- [ ] **Step 7: Commit any fixes and tag**

```bash
git add -A && git commit -m "chore: smoke test fixes"
git tag v0.1.0 -m "Initial capsule scaffold — ready for data asset"
```

---

## Implementation Notes

### Execution order

Tasks 1-6 can be done sequentially (each ~5-15 min). Task 7 (sync script) is the main development effort. Task 8 (dependency audit) can run in parallel with Task 9. Task 11 (path adaptation) is the most labor-intensive — budget 30-60 min for all 8 sections.

### What NOT to do

- Do not modify analysis logic in any script — only paths and config imports
- Do not add features, error handling, or improvements to the analysis code
- Do not create the data asset yet — that's a separate step after CodeOcean quota is resolved
- Do not push to any remote — capsule repo is local only until publication

### Key files to reference during implementation

- Spec: `docs/superpowers/specs/2026-03-22-codeocean-capsule-design.md`
- CITEgeist reference: `../../CITEgeistNeilAnalysis/CITEgeist_CodeOcean/`
- Main repo CLAUDE.md: `CLAUDE.md` (coding standards, compatibility issues)
- Per-section scripts: `analysis/{section}/` in main repo
