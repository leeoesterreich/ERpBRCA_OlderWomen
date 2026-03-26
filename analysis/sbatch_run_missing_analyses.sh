#!/bin/bash
# Master launcher for missing analysis pipelines
# Submits 01, 03, 04 analysis jobs in parallel, then starts polling watcher

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "Missing Analysis Pipeline Launcher"
echo "=============================================="
echo "Directory: $SCRIPT_DIR"
echo "Started: $(date)"
echo ""

# Clean up old markers
MARKER_DIR="$SCRIPT_DIR/.pipeline_markers"
rm -rf "$MARKER_DIR"
mkdir -p "$MARKER_DIR"
echo "Cleaned marker directory: $MARKER_DIR"
echo ""

# Pre-flight checks
echo "=== Pre-flight Checks ==="
REQUIRED_SCRIPTS=(
    "01_human_bulk_rnaseq/run_analysis.sbatch"
    "03_rat_wes/run_analysis.sbatch"
    "04_human_scrnaseq/run_analysis.sbatch"
    "poll_and_validate.sbatch"
)
for script in "${REQUIRED_SCRIPTS[@]}"; do
    if [ ! -f "$script" ]; then
        echo "ERROR: Required script not found: $script"
        exit 1
    fi
    echo "  $script"
done
echo ""

# Submit analysis jobs
echo "=== Submitting Analysis Jobs ==="

JOB1=$(sbatch --parsable 01_human_bulk_rnaseq/run_analysis.sbatch) || { echo "ERROR: Failed to submit 01_human_bulk_rnaseq"; exit 1; }
[ -z "$JOB1" ] && { echo "ERROR: sbatch returned empty job ID for 01_human_bulk_rnaseq"; exit 1; }
echo "01_human_bulk_rnaseq: Job $JOB1"

JOB3=$(sbatch --parsable 03_rat_wes/run_analysis.sbatch) || { echo "ERROR: Failed to submit 03_rat_wes"; exit 1; }
[ -z "$JOB3" ] && { echo "ERROR: sbatch returned empty job ID for 03_rat_wes"; exit 1; }
echo "03_rat_wes: Job $JOB3"

JOB4=$(sbatch --parsable 04_human_scrnaseq/run_analysis.sbatch) || { echo "ERROR: Failed to submit 04_human_scrnaseq"; exit 1; }
[ -z "$JOB4" ] && { echo "ERROR: sbatch returned empty job ID for 04_human_scrnaseq"; exit 1; }
echo "04_human_scrnaseq: Job $JOB4"

echo ""
echo "=== Submitting Polling Watcher ==="

# Submit polling job
POLL_JOB=$(sbatch --parsable poll_and_validate.sbatch) || { echo "ERROR: Failed to submit poll_and_validate"; exit 1; }
[ -z "$POLL_JOB" ] && { echo "ERROR: sbatch returned empty job ID for poll_and_validate"; exit 1; }
echo "poll_and_validate: Job $POLL_JOB"

echo ""
echo "=============================================="
echo "All jobs submitted successfully"
echo "=============================================="
echo "Monitor with: squeue -u $USER"
echo ""
echo "Expected markers:"
echo "  - $MARKER_DIR/01_human_bulk_rnaseq.complete"
echo "  - $MARKER_DIR/03_rat_wes.complete"
echo "  - $MARKER_DIR/04_human_scrnaseq.complete"
echo ""
echo "Validation will run automatically when all complete."
