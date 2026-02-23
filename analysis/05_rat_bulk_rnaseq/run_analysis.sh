#!/bin/bash
# Run rat bulk RNA-seq analysis pipeline
# Usage: bash run_analysis.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p logs outputs/trimmed outputs/aligned outputs/counts

echo "=== Rat Bulk RNA-seq Pipeline ==="
echo "Directory: $SCRIPT_DIR"
echo ""

# Step 0: Trimming
echo "Step 0: Submitting trimming job..."
JOB0=$(sbatch --parsable 00_fastp_trim.sh)
echo "  Job ID: $JOB0"

# Step 1: Alignment (depends on trimming)
echo "Step 1: Submitting alignment job..."
JOB1=$(sbatch --parsable --dependency=afterok:$JOB0 01_alignment.sh)
echo "  Job ID: $JOB1"

# Step 2: HTSeq counting (depends on alignment)
echo "Step 2: Submitting HTSeq job..."
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 02_htseq_count.sh)
echo "  Job ID: $JOB2"

# Step 3-5: R analysis (depends on counting)
echo "Step 3-5: Submitting R analysis job..."
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 run_analysis.sbatch)
echo "  Job ID: $JOB3"

echo ""
echo "=== Jobs submitted ==="
echo "Monitor with: squeue -u \$USER"
echo "Chain: $JOB0 -> $JOB1 -> $JOB2 -> $JOB3"
