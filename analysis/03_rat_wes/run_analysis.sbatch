#!/bin/bash
#SBATCH --job-name=rat_wes
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/wes_%j.out
#SBATCH --error=logs/wes_%j.err

# Rat WES Analysis Pipeline
# Generates COSMIC signatures and oncoplot data

# Load conda environment before strict mode
source /ihome/alee/alc376/.bashrc
conda activate aging_wes || { echo "ERROR: Conda environment 'aging_wes' not found"; exit 1; }

# Enable strict mode after bashrc
set -eo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
echo "Working directory: $SCRIPT_DIR"

mkdir -p logs outputs

echo "=== Step 1: Parse VEP Output ==="
python3 01_parse_vep.py

echo "=== Step 2: COSMIC Signatures ==="
python3 02_cosmic_signatures.py

echo "=== Step 3: Generate Oncoplot Data ==="
python3 03_generate_oncoplot.py

# Write completion marker
MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
echo "Job ${SLURM_JOB_ID:-local} completed at $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$MARKER_DIR/03_rat_wes.complete"
echo "Marker written: $MARKER_DIR/03_rat_wes.complete"

echo "=== Rat WES Analysis Complete ==="
