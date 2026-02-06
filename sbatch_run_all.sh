#!/bin/bash
#SBATCH --job-name=ERpBRCA_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ERpBRCA_OlderWomen Full Analysis Pipeline
# Runs all analyses with biostatistical corrections

set -euo pipefail

# Setup
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
mkdir -p logs

echo "=============================================="
echo "  ERpBRCA OlderWomen Analysis Pipeline"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Load conda
source /ihome/alee/alc376/.bashrc
module load gcc/8.2.0

# Check if environment exists, create if not
if ! conda env list | grep -q "erp_brca_aging"; then
    echo "Creating conda environment..."
    conda env create -f environment.yml
fi

# Activate environment
conda activate erp_brca_aging

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "R version: $(R --version | head -1)"
echo ""

# Run the master script
# Skip download for now since GEO data may not be public yet
# Use --skip-download if data is already present or not yet on GEO
bash scripts/run_all.sh --skip-download

echo ""
echo "=============================================="
echo "  Pipeline Complete"
echo "=============================================="
echo "Finished: $(date)"
