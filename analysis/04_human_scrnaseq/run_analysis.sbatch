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
# Dispatches Stage 1 (preprocessing) and Stage 2 (downstream)
# with SLURM dependency-based checkpointing.
# For independent execution, use run_preprocessing.sbatch or
# run_downstream.sbatch directly.

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq"
cd "$SCRIPT_DIR"
echo "Working directory: $SCRIPT_DIR"

MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"

# --- Stage 1: Preprocessing ---
PREPROC_MARKER="$MARKER_DIR/04_human_scrnaseq_preprocessing.complete"
if [ -f "$PREPROC_MARKER" ]; then
    echo "=== Stage 1 already completed ($(cat "$PREPROC_MARKER")), skipping ==="
    PREPROC_JOB_ID=""
else
    echo "=== Submitting Stage 1: Preprocessing ==="
    PREPROC_JOB_ID=$(sbatch "$SCRIPT_DIR/run_preprocessing.sbatch" | awk '{print $4}')
    echo "Stage 1 submitted as job $PREPROC_JOB_ID"
fi

# --- Stage 2: Downstream ---
DOWNSTREAM_MARKER="$MARKER_DIR/04_human_scrnaseq.complete"
if [ -f "$DOWNSTREAM_MARKER" ]; then
    echo "=== Stage 2 already completed ($(cat "$DOWNSTREAM_MARKER")), skipping ==="
else
    SBATCH_ARGS=""
    if [ -n "$PREPROC_JOB_ID" ]; then
        SBATCH_ARGS="--dependency=afterok:$PREPROC_JOB_ID"
        echo "Stage 2 will wait for Stage 1 (job $PREPROC_JOB_ID) to complete"
    fi
    echo "=== Submitting Stage 2: Downstream ==="
    DOWNSTREAM_JOB_ID=$(sbatch $SBATCH_ARGS "$SCRIPT_DIR/run_downstream.sbatch" | awk '{print $4}')
    echo "Stage 2 submitted as job $DOWNSTREAM_JOB_ID"
fi

echo "=== Pipeline dispatch complete ==="
