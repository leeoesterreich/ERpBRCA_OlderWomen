#!/bin/bash
#SBATCH --job-name=fig7_r_only
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fig7_r_only_%j.out
#SBATCH --error=logs/fig7_r_only_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Figure 7 R scripts only (CellPhoneDB already completed)
set -eo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname $(readlink -f $0))}"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "=============================================="
echo "Figure 7 R Analysis"
echo "Job ID: $SLURM_JOB_ID"
echo "Start: $(date)"
echo "=============================================="

source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

# Step 1: Run expanded pathway analysis
echo ""
echo "Step 1: Running expanded pathway analysis..."
Rscript 15_multicelltype_pathway.R

# Step 2: Generate CellPhoneDB dot plot
echo ""
echo "Step 2: Generating CellPhoneDB dot plot..."
Rscript 17_cellphonedb_dotplot.R

echo ""
echo "=============================================="
echo "Completed: $(date)"
echo "=============================================="

# List outputs
echo "Output files:"
ls -la figures/fig7*.png 2>/dev/null || echo "  No figure files found"
