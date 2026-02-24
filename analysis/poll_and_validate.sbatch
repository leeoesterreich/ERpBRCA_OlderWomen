#!/bin/bash
#SBATCH --job-name=poll_validate
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 06:00:00
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/poll_validate_%j.out
#SBATCH --error=logs/poll_validate_%j.err

# Polling watcher that triggers validation when all analyses complete

set -eo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

MARKER_DIR="$SCRIPT_DIR/.pipeline_markers"
POLL_INTERVAL=60  # seconds
MAX_WAIT=21600    # 6 hours in seconds

REQUIRED_MARKERS=(
    "01_human_bulk_rnaseq.complete"
    "03_rat_wes.complete"
    "04_human_scrnaseq.complete"
)

echo "=============================================="
echo "Polling Watcher for Analysis Completion"
echo "=============================================="
echo "Marker directory: $MARKER_DIR"
echo "Poll interval: ${POLL_INTERVAL}s"
echo "Max wait: $((MAX_WAIT / 3600))h"
echo "Required markers: ${#REQUIRED_MARKERS[@]}"
echo ""

wait_start=$(date +%s)

while true; do
    elapsed=$(($(date +%s) - wait_start))

    if [ $elapsed -ge $MAX_WAIT ]; then
        echo "ERROR: Timeout after $((elapsed / 3600))h waiting for markers"
        echo "Missing markers:"
        for marker in "${REQUIRED_MARKERS[@]}"; do
            if [ ! -f "$MARKER_DIR/$marker" ]; then
                echo "  - $marker"
            fi
        done
        exit 1
    fi

    # Check all markers
    all_complete=true
    for marker in "${REQUIRED_MARKERS[@]}"; do
        if [ ! -f "$MARKER_DIR/$marker" ]; then
            all_complete=false
            break
        fi
    done

    if [ "$all_complete" = true ]; then
        echo ""
        echo "All markers present! Starting validation..."
        echo "Elapsed time: $((elapsed / 60)) minutes"
        break
    fi

    # Status update every 5 minutes
    if [ $((elapsed % 300)) -lt $POLL_INTERVAL ]; then
        echo "[$(date '+%H:%M:%S')] Waiting... ($((elapsed / 60))m elapsed)"
        for marker in "${REQUIRED_MARKERS[@]}"; do
            if [ -f "$MARKER_DIR/$marker" ]; then
                echo "  ✓ $marker"
            else
                echo "  ○ $marker (pending)"
            fi
        done
    fi

    sleep $POLL_INTERVAL
done

echo ""
echo "=============================================="
echo "Running Validation Pipeline"
echo "=============================================="

# Activate environment for validation
source /ihome/alee/alc376/.bashrc
conda activate erp_brca_aging || { echo "ERROR: Failed to activate conda environment"; exit 1; }

# Change to project root and run validation
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

bash scripts/validation/run_validation.sh

echo ""
echo "=============================================="
echo "Validation Complete"
echo "=============================================="
echo "Report: $PROJECT_ROOT/validation/report/validation_report.md"
