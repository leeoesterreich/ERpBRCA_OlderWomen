# Spatial Biopsies Section Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `analysis/06_spatial_biopsies/` section to ERpBRCA_OlderWomen by copying and adapting CITEgeist scripts 12 and 13.

**Architecture:** Minimal copy of two Python scripts plus extracted utility and figure config into a new analysis section, with symlinked data and updated sbatch/env references.

**Tech Stack:** Python (scanpy, gseapy, seaborn, psutil), SLURM

**Spec:** `docs/superpowers/specs/2026-03-15-spatial-biopsies-design.md`

---

## Chunk 1: Directory Structure and Dependencies

### Task 1: Create directory structure and data symlink

**Files:**
- Create: `analysis/06_spatial_biopsies/data/` (directory)
- Create: `analysis/06_spatial_biopsies/logs/` (directory)
- Create: `analysis/06_spatial_biopsies/figures/immune_secretion/` (directory)
- Create: `analysis/06_spatial_biopsies/figures/immune_pathways/enrichment/` (directory)
- Create: `analysis/06_spatial_biopsies/figures/immune_pathways/prerank/` (directory)
- Create: `analysis/06_spatial_biopsies/figures/immune_pathways/dotplots/` (directory)
- Create: `analysis/06_spatial_biopsies/data/biopsy_adatas.pkl` (symlink)

- [ ] **Step 1: Create all directories**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
mkdir -p analysis/06_spatial_biopsies/{data,logs}
mkdir -p analysis/06_spatial_biopsies/figures/{immune_secretion,immune_pathways/{enrichment,prerank,dotplots}}
```

- [ ] **Step 2: Create data symlink**

```bash
ln -s /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/data/biopsy_adatas.pkl \
  analysis/06_spatial_biopsies/data/biopsy_adatas.pkl
```

- [ ] **Step 3: Verify symlink resolves**

```bash
ls -la analysis/06_spatial_biopsies/data/biopsy_adatas.pkl
```
Expected: symlink pointing to CITEgeist data file.

- [ ] **Step 4: Add .gitkeep files to empty directories**

```bash
touch analysis/06_spatial_biopsies/logs/.gitkeep
touch analysis/06_spatial_biopsies/figures/immune_secretion/.gitkeep
touch analysis/06_spatial_biopsies/figures/immune_pathways/enrichment/.gitkeep
touch analysis/06_spatial_biopsies/figures/immune_pathways/prerank/.gitkeep
touch analysis/06_spatial_biopsies/figures/immune_pathways/dotplots/.gitkeep
```

- [ ] **Step 5: Commit directory structure**

```bash
git add analysis/06_spatial_biopsies/
git commit -m "chore: create 06_spatial_biopsies directory structure with data symlink"
```

### Task 2: Create utils.py (extracted)

**Files:**
- Create: `analysis/06_spatial_biopsies/utils.py`

- [ ] **Step 1: Write utils.py with only log_memory_usage**

```python
"""Utility functions for spatial biopsy analysis."""

import os
import logging
import psutil


def log_memory_usage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1024 / 1024
    logging.info(f"Memory usage: {mem:.2f} MB")
```

- [ ] **Step 2: Commit**

```bash
git add analysis/06_spatial_biopsies/utils.py
git commit -m "feat: add utils.py with log_memory_usage for spatial biopsies"
```

### Task 3: Copy and update figure_config.py

**Files:**
- Create: `analysis/06_spatial_biopsies/figure_config.py`

Source: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/figure_config.py`

- [ ] **Step 1: Copy figure_config.py**

```bash
cp /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/figure_config.py \
  analysis/06_spatial_biopsies/figure_config.py
```

- [ ] **Step 2: Update docstring**

Change line 1 from:
```python
"""Centralized figure configuration for CITEgeist analysis scripts.
```
To:
```python
"""Centralized figure configuration for spatial biopsy analysis scripts.
```

- [ ] **Step 3: Commit**

```bash
git add analysis/06_spatial_biopsies/figure_config.py
git commit -m "feat: add figure_config.py for spatial biopsies"
```

### Task 4: Verify Python dependencies in erp_brca_aging

- [ ] **Step 1: Check which packages are already installed**

```bash
conda run -p /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging \
  python -c "import scanpy; import gseapy; import seaborn; import psutil; print('All dependencies available')" 2>&1
```

- [ ] **Step 2: Install any missing packages**

If any imports fail, install them:

```bash
conda run -p /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging \
  pip install <missing_package>
```

No commit needed — env changes are not tracked in git.

---

## Chunk 2: Copy and Adapt Analysis Scripts

### Task 5: Create 01_immune_secretion.py

**Files:**
- Create: `analysis/06_spatial_biopsies/01_immune_secretion.py`

Source: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/12_immune_secretion.py`

- [ ] **Step 1: Copy source script**

```bash
cp /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/12_immune_secretion.py \
  analysis/06_spatial_biopsies/01_immune_secretion.py
```

- [ ] **Step 2: Update log filename**

Change line 21 from:
```python
        logging.FileHandler('logs/12_immune_secretion.log'),
```
To:
```python
        logging.FileHandler('logs/01_immune_secretion.log'),
```

- [ ] **Step 3: Verify no other changes needed**

The script uses only relative paths (`data/`, `logs/`, `figures/`) which are correct
for the new location. The `from utils import log_memory_usage` import will resolve
since `utils.py` is in the same directory and sbatch will `cd` there.

- [ ] **Step 4: Commit**

```bash
git add analysis/06_spatial_biopsies/01_immune_secretion.py
git commit -m "feat: add 01_immune_secretion.py for spatial biopsies"
```

### Task 6: Create 02_immune_pathways.py

**Files:**
- Create: `analysis/06_spatial_biopsies/02_immune_pathways.py`

Source: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/13_immune_pathways.py`

- [ ] **Step 1: Copy source script**

```bash
cp /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/13_immune_pathways.py \
  analysis/06_spatial_biopsies/02_immune_pathways.py
```

- [ ] **Step 2: Update log filename**

Change line 34 from:
```python
        logging.FileHandler('logs/13_immune_pathways.log'),
```
To:
```python
        logging.FileHandler('logs/02_immune_pathways.log'),
```

- [ ] **Step 3: Update GSEA seed**

Change the `seed=42` in the `gp.prerank()` call (around line 163) from:
```python
                            seed=42,
```
To:
```python
                            seed=12345,
```

- [ ] **Step 4: Verify no other changes needed**

Same relative path logic as script 01 — all paths are relative and correct.

- [ ] **Step 5: Commit**

```bash
git add analysis/06_spatial_biopsies/02_immune_pathways.py
git commit -m "feat: add 02_immune_pathways.py for spatial biopsies"
```

---

## Chunk 3: Sbatch Runner and Project Updates

### Task 7: Create run_analysis.sbatch

**Files:**
- Create: `analysis/06_spatial_biopsies/run_analysis.sbatch`

Follow the pattern from `analysis/01_human_bulk_rnaseq/run_analysis.sbatch`.

- [ ] **Step 1: Write run_analysis.sbatch**

```bash
#!/bin/bash
#SBATCH --job-name=spatial_biopsies
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/06_spatial_biopsies/logs/spatial_biopsies_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/06_spatial_biopsies/logs/spatial_biopsies_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Spatial Biopsies Analysis Pipeline
# Runs immune secretion spatial plots and immune pathway enrichment

set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/06_spatial_biopsies"
cd "$SCRIPT_DIR"

# Create directories
mkdir -p logs figures/immune_secretion
mkdir -p figures/immune_pathways/{enrichment,prerank,dotplots}

echo "=============================================="
echo "Spatial Biopsies Analysis Pipeline"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Working directory: $SCRIPT_DIR"
echo "Start time: $(date)"
echo ""

# Load conda and activate environment
set +u
source ~/.bashrc
eval "$(conda shell.bash hook)"
set -u
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

echo "Python version: $(python --version)"
echo ""

# Verify data symlink
if [ ! -f data/biopsy_adatas.pkl ]; then
    echo "ERROR: data/biopsy_adatas.pkl not found (symlink broken?)"
    exit 1
fi

# Run analysis scripts in order
echo "[1/2] Immune secretion spatial plots..."
python 01_immune_secretion.py 2>&1

echo "[2/2] Immune pathway enrichment..."
python 02_immune_pathways.py "$@" 2>&1

# Write completion marker
MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
touch "$MARKER_DIR/06_spatial_biopsies.complete"
echo "Marker written: $MARKER_DIR/06_spatial_biopsies.complete"

echo ""
echo "=============================================="
echo "Spatial Biopsies Analysis Complete"
echo "=============================================="
echo "End time: $(date)"
echo "Figures: $SCRIPT_DIR/figures/"
```

Note: `"$@"` passes through any extra args (e.g. `--figures-only`) to script 02.

- [ ] **Step 2: Commit**

```bash
git add analysis/06_spatial_biopsies/run_analysis.sbatch
git commit -m "feat: add run_analysis.sbatch for spatial biopsies"
```

### Task 8: Create run_analysis.sh

**Files:**
- Create: `analysis/06_spatial_biopsies/run_analysis.sh`

- [ ] **Step 1: Write run_analysis.sh**

```bash
#!/bin/bash
# analysis/06_spatial_biopsies/run_analysis.sh
# Runs spatial biopsies analysis pipeline
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Spatial Biopsies Analysis Pipeline ==="
echo "Working directory: $SCRIPT_DIR"
echo ""

# Create directories
mkdir -p logs figures/immune_secretion
mkdir -p figures/immune_pathways/{enrichment,prerank,dotplots}

# Verify data
if [ ! -f data/biopsy_adatas.pkl ]; then
    echo "ERROR: data/biopsy_adatas.pkl not found"
    exit 1
fi

echo "[1/2] Immune secretion spatial plots..."
python 01_immune_secretion.py

echo "[2/2] Immune pathway enrichment..."
python 02_immune_pathways.py "$@"

echo ""
echo "=== Spatial Biopsies Analysis Complete ==="
echo "Figures: $SCRIPT_DIR/figures/"
```

- [ ] **Step 2: Make executable**

```bash
chmod +x analysis/06_spatial_biopsies/run_analysis.sh
```

- [ ] **Step 3: Commit**

```bash
git add analysis/06_spatial_biopsies/run_analysis.sh
git commit -m "feat: add run_analysis.sh for spatial biopsies"
```

### Task 9: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add 06_spatial_biopsies to project structure**

In the `## Project Structure` section, add after `05_rat_bulk_rnaseq/`:

```
├── 06_spatial_biopsies/   # Spatial biopsy immune analysis (HCC22-088)
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: add 06_spatial_biopsies to CLAUDE.md project structure"
```
