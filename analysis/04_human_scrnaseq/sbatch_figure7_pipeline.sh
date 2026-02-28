#!/bin/bash
#SBATCH --job-name=fig7_pipeline
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/fig7_pipeline_%j.out
#SBATCH --error=logs/fig7_pipeline_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Figure 7 regeneration pipeline
# 1. Run CellPhoneDB prep (R) - fixes empty cell types
# 2. Run CellPhoneDB analysis (Python)
# 3. Run pathway analysis (R) - with expanded pathways

set -eo pipefail

# Get script directory
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname $(readlink -f $0))}"
cd "$SCRIPT_DIR"
mkdir -p logs

echo "=============================================="
echo "Figure 7 Pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"
echo "Directory: $(pwd)"
echo "=============================================="

# Load conda
source ~/.bashrc
CONDA_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging"

# Step 1: Run CellPhoneDB prep (R)
echo ""
echo "Step 1: Preparing CellPhoneDB input files..."
conda activate "$CONDA_ENV"
Rscript 09_cellphonedb_prep.R

# Step 2: Run CellPhoneDB analysis (Python)
echo ""
echo "Step 2: Running CellPhoneDB analysis..."
# Use base conda for CellPhoneDB (installed there)
conda deactivate
conda activate base
python 16_run_cellphonedb.py

# Step 3: Run pathway analysis (R)
echo ""
echo "Step 3: Running expanded pathway analysis..."
conda activate "$CONDA_ENV"
Rscript 15_multicelltype_pathway.R

# Step 4: Generate CellPhoneDB dot plot (R)
echo ""
echo "Step 4: Generating CellPhoneDB dot plot..."
Rscript 17_cellphonedb_dotplot.R

echo ""
echo "=============================================="
echo "Pipeline completed: $(date)"
echo "=============================================="

# List outputs
echo ""
echo "Output files:"
ls -la figures/fig7*.png 2>/dev/null || echo "  No figure files found"
ls -la outputs/multicelltype_pathway_scores.csv 2>/dev/null || echo "  No pathway scores found"
ls -la outputs/cellphonedb/results_*/pvalues.csv 2>/dev/null || echo "  No CellPhoneDB results found"
