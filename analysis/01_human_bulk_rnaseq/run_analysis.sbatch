#!/bin/bash
#SBATCH --job-name=human_bulk_rnaseq
#SBATCH --output=logs/bulk_rnaseq_%j.out
#SBATCH --error=logs/bulk_rnaseq_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Human Bulk RNA-seq Analysis Pipeline
# Runs preprocessing, GSVA, PROGENy, correlations, and visualization

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Create log directory
mkdir -p logs outputs

echo "=============================================="
echo "Human Bulk RNA-seq Analysis Pipeline"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Working directory: $SCRIPT_DIR"
echo "Start time: $(date)"
echo ""

# Load conda and activate environment
source ~/.bashrc
conda activate erp_bulk_rnaseq || {
    echo "Conda environment not found, attempting to create..."
    conda env create -f ../../envs/bulk_rnaseq.yml
    conda activate erp_bulk_rnaseq
}

echo "R version: $(R --version | head -1)"
echo ""

# Run analysis scripts in order
echo "[1/5] Preprocessing..."
Rscript 01_preprocess.R 2>&1

echo "[2/5] Running GSVA..."
Rscript 02_run_gsva.R 2>&1

echo "[3/5] Running PROGENy..."
Rscript 03_run_progeny.R 2>&1

echo "[4/5] Calculating correlations..."
Rscript 04_correlations.R 2>&1

echo "[5/5] Generating visualizations..."
Rscript 05_visualize.R 2>&1

echo ""
echo "=============================================="
echo "Human Bulk RNA-seq Analysis Complete"
echo "=============================================="
echo "End time: $(date)"
echo "Outputs: $SCRIPT_DIR/outputs/"
