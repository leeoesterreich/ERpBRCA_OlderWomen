# Pipeline Checkpointing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split the monolithic scRNA-seq pipeline into two independently-submittable sbatch scripts with checkpoint-based dependency so partial failures don't require restarting from Step 01.

**Architecture:** Create `run_preprocessing.sbatch` (Steps 00b-05, 256G/12h) and `run_downstream.sbatch` (Steps 06-11, 16G/4h). Replace `run_analysis.sbatch` with an orchestrator that checks completion markers before dispatching each stage. R scripts remain untouched.

**Tech Stack:** Bash, SLURM sbatch, R scripts (unchanged)

## Global Constraints

- HTC cluster: `#SBATCH --cluster=htc` and `#SBATCH --partition=htc` — never use `smp` or `gpu`
- Conda on HTC: `export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin:$PATH"` then `source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh`
- Environment: `conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging`
- Mail on failure: `#SBATCH --mail-type=FAIL` and `#SBATCH --mail-user=alc376@pitt.edu`
- Working directory: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq`
- Marker directory: `analysis/.pipeline_markers/` (one level up from `04_human_scrnaseq/`)
- Don't modify any R scripts — only sbatch files change
- `set -eo pipefail` in all scripts (after conda activation)

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `run_preprocessing.sbatch` | CREATE | Stage 1: Steps 00b-05, writes preprocessing marker |
| `run_downstream.sbatch` | CREATE | Stage 2: Steps 06-11, writes final completion marker |
| `run_analysis.sbatch` | MODIFY | Orchestrator: checks markers, submits stages |

---

### Task 1: Create run_preprocessing.sbatch

**Files:**
- Create: `analysis/04_human_scrnaseq/run_preprocessing.sbatch`
- Reference: `analysis/04_human_scrnaseq/run_analysis.sbatch` (lines 1-54)

**Interfaces:**
- Consumes: `../../data/human_scrnaseq/SeuratObj_Xu2024_HRpos_AfterQCSCT.rds` (37GB shared file)
- Produces: `outputs/seurat_annotated.rds` (640MB checkpoint), `.pipeline_markers/04_human_scrnaseq_preprocessing.complete`

- [ ] **Step 1: Write run_preprocessing.sbatch**

```bash
#!/bin/bash
#SBATCH --job-name=scrnaseq_preproc
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/preprocessing_%j.out
#SBATCH --error=logs/preprocessing_%j.err

# Stage 1: Load, annotate, cell fractions, gene expression
# Produces: outputs/seurat_annotated.rds (checkpoint for Stage 2)

# Load conda environment on HTC cluster
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin:$PATH"
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

# Now enable strict mode
set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq"
cd "$SCRIPT_DIR"
echo "Working directory: $SCRIPT_DIR"

mkdir -p logs outputs outputs/preprocessing figures

# Step 0b: Preprocess raw data (skip if already done)
PREPROCESS_FILE="../../data/human_scrnaseq/SeuratObj_Xu2024_HRpos_AfterQCSCT.rds"
if [ -f "$PREPROCESS_FILE" ]; then
    echo "=== Step 0b: Preprocessed data exists, skipping ==="
else
    echo "=== Step 0b: Preprocess raw data ==="
    Rscript 00b_preprocess_seurat.R
fi

echo "=== Step 1: Load and subset data ==="
Rscript 01_load_subset_data.R

# Note: Harmony integration removed (not used per methods)

echo "=== Step 3: Cell type annotation ==="
Rscript 03_cell_type_annotation.R

echo "=== Step 4: Cell fractions ==="
Rscript 04_cell_fractions.R

echo "=== Step 5: Gene expression ==="
Rscript 05_gene_expression_violin.R

# Validate checkpoint file was created
if [ ! -f "outputs/seurat_annotated.rds" ]; then
    echo "ERROR: seurat_annotated.rds not created — Stage 1 failed"
    exit 1
fi

# Write preprocessing completion marker
MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
echo "Stage 1 (preprocessing) Job ${SLURM_JOB_ID:-local} completed at $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$MARKER_DIR/04_human_scrnaseq_preprocessing.complete"
echo "Marker written: $MARKER_DIR/04_human_scrnaseq_preprocessing.complete"

echo "=== Preprocessing complete ==="
```

- [ ] **Step 2: Verify file content**

```bash
cat analysis/04_human_scrnaseq/run_preprocessing.sbatch | head -5
```
Expected: Shebang + first 4 SBATCH directives visible.

- [ ] **Step 3: Set executable permission**

```bash
chmod +x analysis/04_human_scrnaseq/run_preprocessing.sbatch
```

- [ ] **Step 4: Commit**

```bash
git add analysis/04_human_scrnaseq/run_preprocessing.sbatch
git commit -m "feat(pipeline): add run_preprocessing.sbatch for Stage 1 (Steps 00b-05)"
```

---

### Task 2: Create run_downstream.sbatch

**Files:**
- Create: `analysis/04_human_scrnaseq/run_downstream.sbatch`
- Reference: `analysis/04_human_scrnaseq/run_analysis.sbatch` (lines 55-141)

**Interfaces:**
- Consumes: `outputs/seurat_annotated.rds` (from Stage 1), `../../data/human_scrnaseq/SeuratObj_Xu2024_HRpos_AfterQCSCT.rds`
- Produces: GSVA figures, PROGENy/WCSEA/CellPhoneDB outputs, `.pipeline_markers/04_human_scrnaseq.complete`

- [ ] **Step 1: Write run_downstream.sbatch**

Copy the conda activation block (lines 18-28) from `run_analysis.sbatch`, then Steps 06-11 (lines 55-141). Key differences from the original:

- `--mem=16G` and `--cpus-per-task=4` (downstream doesn't need 256G)
- `--time 4:00:00` (downstream is fast)
- Input validation: check `seurat_annotated.rds` exists before starting
- Job name: `scrnaseq_downstream`
- Log pattern: `logs/downstream_%j.out`
- Uses the existing completion marker path (`04_human_scrnaseq.complete`) for backwards compatibility

```bash
#!/bin/bash
#SBATCH --job-name=scrnaseq_downstream
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/downstream_%j.out
#SBATCH --error=logs/downstream_%j.err

# Stage 2: GSVA, PROGENy, WCSEA, CellPhoneDB prep, TAM export
# Depends on: outputs/seurat_annotated.rds from Stage 1

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin:$PATH"
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq"
cd "$SCRIPT_DIR"
echo "Working directory: $SCRIPT_DIR"

mkdir -p logs outputs outputs/preprocessing figures

# Validate Stage 1 checkpoint exists
if [ ! -f "outputs/seurat_annotated.rds" ]; then
    echo "ERROR: outputs/seurat_annotated.rds not found. Run Stage 1 (run_preprocessing.sbatch) first."
    exit 1
fi

echo "=== Step 6: GSVA (pseudobulk only) ==="
Rscript 06_run_gsva.R --mode=pseudobulk

echo "=== Step 7: PROGENy ==="
Rscript 07_run_progeny.R

echo "=== Step 8: WCSEA ==="
Rscript 08_run_wcsea.R

echo "=== Step 9: CellPhoneDB prep ==="
Rscript 09_cellphonedb_prep.R

echo "=== Step 10: Export TAM analysis CSV ==="
Rscript -e '
tryCatch({
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(dplyr))

  if (!file.exists("outputs/seurat_annotated.rds")) {
    stop("seurat_annotated.rds not found")
  }

  seurat <- readRDS("outputs/seurat_annotated.rds")

  if (!"CellTypeAnnot" %in% colnames(seurat@meta.data)) {
    stop("Column CellTypeAnnot not found in Seurat metadata")
  }

  mac_cells <- WhichCells(seurat, expression = CellTypeAnnot == "Macrophage")

  if (length(mac_cells) == 0) {
    warning("No macrophage cells found, creating empty TAM analysis")
    tam_metrics <- data.frame(sample = character(), age_group = character(), n_cells = integer())
  } else {
    mac_subset <- subset(seurat, cells = mac_cells)

    tam_metrics <- mac_subset@meta.data %>%
      group_by(orig.ident) %>%
      summarise(
        n_cells = n(),
        .groups = "drop"
      )

    if ("AgeGroup" %in% colnames(mac_subset@meta.data)) {
      age_mapping <- mac_subset@meta.data %>%
        select(orig.ident, AgeGroup) %>%
        distinct()
      tam_metrics <- left_join(tam_metrics, age_mapping, by = "orig.ident")
    }
  }

  write.csv(tam_metrics, "outputs/tam_analysis.csv", row.names = FALSE)
  cat("TAM analysis exported:", nrow(tam_metrics), "samples\n")
}, error = function(e) {
  message(paste("ERROR in TAM export:", e$message))
  quit(status = 1)
})
' || { echo "ERROR: TAM export failed"; exit 1; }

echo "=== Step 11: Ensure PNG figures exist ==="
for pdf in outputs/*.pdf; do
    [ -f "$pdf" ] || continue
    base=$(basename "$pdf" .pdf)
    png="figures/${base}.png"
    if [ ! -f "$png" ]; then
        echo "Converting $pdf -> $png"
        env -i PATH=/usr/bin:/bin pdftoppm -png -r 300 -singlefile "$pdf" "figures/${base}" 2>/dev/null \
            || env -i PATH=/usr/bin:/bin convert -density 300 "$pdf" "$png" 2>/dev/null \
            || echo "WARNING: Could not convert $pdf"
    fi
done

MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
echo "Stage 2 (downstream) Job ${SLURM_JOB_ID:-local} completed at $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$MARKER_DIR/04_human_scrnaseq.complete"
echo "Marker written: $MARKER_DIR/04_human_scrnaseq.complete"

echo "=== Downstream analysis complete ==="
```

- [ ] **Step 2: Set executable permission**

```bash
chmod +x analysis/04_human_scrnaseq/run_downstream.sbatch
```

- [ ] **Step 3: Commit**

```bash
git add analysis/04_human_scrnaseq/run_downstream.sbatch
git commit -m "feat(pipeline): add run_downstream.sbatch for Stage 2 (Steps 06-11)"
```

---

### Task 3: Rewrite run_analysis.sbatch as Orchestrator

**Files:**
- Modify: `analysis/04_human_scrnaseq/run_analysis.sbatch`

**Interfaces:**
- Consumes: `run_preprocessing.sbatch`, `run_downstream.sbatch`
- Produces: Same completion marker as before (backwards compatible)

- [ ] **Step 1: Replace run_analysis.sbatch with orchestrator**

```bash
#!/bin/bash
#SBATCH --job-name=human_scrnaseq
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/analysis_%j.out
#SBATCH --error=logs/analysis_%j.err

# Human scRNA-seq analysis pipeline orchestrator
# Dispatches Stage 1 (preprocessing) and Stage 2 (downstream) with checkpointing

set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq"
cd "$SCRIPT_DIR"
echo "Working directory: $SCRIPT_DIR"

MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"

# --- Stage 1: Preprocessing ---
PREPROC_MARKER="$MARKER_DIR/04_human_scrnaseq_preprocessing.complete"
if [ -f "$PREPROC_MARKER" ]; then
    echo "=== Stage 1 already completed ($(cat "$PREPROC_MARKER")), skipping ==="
else
    echo "=== Running Stage 1: Preprocessing ==="
    bash "$SCRIPT_DIR/run_preprocessing.sbatch" &
    PREPROC_PID=$!
    wait $PREPROC_PID || { echo "ERROR: Stage 1 failed"; exit 1; }
    echo "=== Stage 1 complete ==="
fi

# --- Stage 2: Downstream ---
DOWNSTREAM_MARKER="$MARKER_DIR/04_human_scrnaseq.complete"
if [ -f "$DOWNSTREAM_MARKER" ]; then
    echo "=== Stage 2 already completed ($(cat "$DOWNSTREAM_MARKER")), skipping ==="
else
    echo "=== Running Stage 2: Downstream ==="
    bash "$SCRIPT_DIR/run_downstream.sbatch" &
    DOWNSTREAM_PID=$!
    wait $DOWNSTREAM_PID || { echo "ERROR: Stage 2 failed"; exit 1; }
    echo "=== Stage 2 complete ==="
fi

echo "=== All stages complete ==="
```

**Note:** This orchestrator runs the stage scripts locally (not via sbatch) because we're already on a compute node. The SBATCH directives in the stage scripts are ignored by `bash`. The conda activation inside each stage script still works correctly.

- [ ] **Step 2: Commit**

```bash
git add analysis/04_human_scrnaseq/run_analysis.sbatch
git commit -m "refactor(pipeline): convert run_analysis.sbatch to checkpoint-based orchestrator"
```

---

### Task 4: Submit and Test Stage 2 Standalone

**Files:**
- No file changes — testing only

**Interfaces:**
- Consumes: `outputs/seurat_annotated.rds` (already exists from the previous run)
- Produces: GSVA outputs, completion marker

- [ ] **Step 1: Push changes to remote**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
git push origin main
```

- [ ] **Step 2: Submit run_downstream.sbatch to HTC cluster**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq
sbatch run_downstream.sbatch
```
Expected: SLURM job ID returned, e.g. `Submitted batch job 10338430`

- [ ] **Step 3: Monitor job**

Use the `/babysit` skill to monitor the submitted job.

- [ ] **Step 4: Verify outputs**

When job completes, check:
```bash
ls -la figures/gsva_heatmap_rotated_pseudobulk.*
ls -la figures/gsva_per_patient_comparison_pseudobulk.*
ls -la outputs/gsva_comparison_stats_pseudobulk.csv
cat ../.pipeline_markers/04_human_scrnaseq.complete
```

Success criteria: All three output paths exist and have non-zero size. Completion marker exists.

- [ ] **Step 5: Commit any new outputs**

```bash
git add figures/gsva_* figures/gsva_* outputs/gsva_*
git commit -m "results(gsva): add rotated HSD17B7 +/ heatmap and per-patient comparison"
```

---

## Self-Review

**1. Spec coverage:**
- Two checkpointed sbatch scripts? ✓ (Tasks 1, 2)
- Orchestrator with marker checking? ✓ (Task 3)
- Backwards-compatible completion marker? ✓ (`04_human_scrnaseq.complete` unchanged)
- Resource differentiation (256G vs 16G)? ✓
- R scripts untouched? ✓
- Testing strategy (Stage 2 standalone)? ✓ (Task 4)

**2. Placeholder scan:** No TBDs, no "add validation," no vague steps. All code blocks are complete.

**3. Type consistency:** Marker file paths consistent across all three scripts. Conda paths match. SBATCH directives match spec. No signature mismatches (all bash, no function interfaces between tasks).

**One risk to flag:** The orchestrator (Task 3) runs stage scripts via `bash` (not `sbatch`) since it's already on a compute node. This means the SBATCH directives in stage scripts are ignored — the orchestrator's own SBATCH directives govern resources. This is correct behavior (the orchestrator requests 256G/12h, which is enough for both stages). If Stage 2 starts failing with memory errors, we may need to re-submit Stage 2 independently with sufficient memory — that's the whole point of the split.
