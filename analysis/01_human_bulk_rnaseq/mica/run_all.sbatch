#!/bin/bash
#SBATCH --job-name=mica_analysis
#SBATCH --partition=htc
#SBATCH --time=4:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=10
#SBATCH --output=logs/mica_%j.out
#SBATCH --error=logs/mica_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

# Activate environment
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
cd "$SCRIPT_DIR"

echo "=== MICA Analysis Pipeline ==="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"

# Run pipeline sequentially
echo "Step 1: Preparing data..."
Rscript 01_prep_data.R

echo "Step 2: Running GSVA..."
Rscript 02_run_gsva.R

echo "Step 3: Running MICA..."
Rscript 03_run_mica.R

echo "Step 4: Generating Figure 4..."
Rscript 04_fig4_estrogen.R

echo "Step 5: Generating Figure 5..."
Rscript 05_fig5_immune.R

echo "=== Pipeline complete ==="
echo "End time: $(date)"
