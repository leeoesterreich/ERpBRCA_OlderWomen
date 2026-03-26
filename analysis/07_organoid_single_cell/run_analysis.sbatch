#!/bin/bash
#SBATCH --job-name=organoid_sc
#SBATCH --output=logs/organoid_sc_%j.out
#SBATCH --error=logs/organoid_sc_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Organoid Single-Cell Analysis Pipeline
# Runs QC through proliferation analysis (scripts 02-11)

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "Organoid Single-Cell Analysis Pipeline"
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
conda activate erp_brca_aging

echo "Python version: $(python --version)"
echo ""

# Run full pipeline
bash run_analysis.sh

# Write completion marker (matches section 06 convention — marker in sbatch, not sh)
MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
touch "$MARKER_DIR/07_organoid_single_cell.complete"
echo "Marker written: $MARKER_DIR/07_organoid_single_cell.complete"

echo ""
echo "=============================================="
echo "Organoid Single-Cell Analysis Complete"
echo "=============================================="
echo "End time: $(date)"
