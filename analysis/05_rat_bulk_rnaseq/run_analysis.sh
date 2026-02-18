#!/bin/bash
# Run rat bulk RNA-seq analysis pipeline
# Usage: bash run_analysis.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Rat Bulk RNA-seq Pipeline ==="
echo "Directory: $SCRIPT_DIR"

# Step 1: Alignment (submit as SLURM job)
echo "Step 1: Submitting alignment job..."
JOB1=$(sbatch --parsable 01_alignment.sh)
echo "  Job ID: $JOB1"

# Step 2: HTSeq counting (depends on alignment)
echo "Step 2: Submitting HTSeq job..."
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 02_htseq_count.sh)
echo "  Job ID: $JOB2"

# Step 3-4: R analysis (submit after counting)
echo "Step 3-4: Submitting R analysis job..."
sbatch --dependency=afterok:$JOB2 run_analysis.sbatch

echo "=== Jobs submitted ==="
