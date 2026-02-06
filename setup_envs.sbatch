#!/bin/bash
#SBATCH --job-name=setup_envs
#SBATCH --output=logs/setup_envs_%j.out
#SBATCH --error=logs/setup_envs_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Setup conda environments for ERp BRCA analyses

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p logs

echo "=============================================="
echo "Setting up conda environments"
echo "=============================================="
echo "Start time: $(date)"
echo ""

source ~/.bashrc

# Prefer mamba for faster solving
CONDA_CMD="mamba"
if ! command -v mamba &>/dev/null; then
    CONDA_CMD="conda"
fi
echo "Using: $CONDA_CMD"

# Create bulk RNA-seq environment
if conda env list | grep -q "erp_bulk_rnaseq"; then
    echo "erp_bulk_rnaseq environment already exists"
else
    echo "Creating erp_bulk_rnaseq environment..."
    $CONDA_CMD env create -f envs/bulk_rnaseq.yml
    echo "erp_bulk_rnaseq created successfully"
fi

echo ""

# Create snRNA-seq environment
if conda env list | grep -q "erp_snrnaseq"; then
    echo "erp_snrnaseq environment already exists"
else
    echo "Creating erp_snrnaseq environment..."
    $CONDA_CMD env create -f envs/snrnaseq.yml
    echo "erp_snrnaseq created successfully"
fi

echo ""
echo "=============================================="
echo "Environment setup complete"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "To run the analyses, submit:"
echo "  sbatch analysis/01_human_bulk_rnaseq/run_analysis.sbatch"
echo "  sbatch analysis/02_rat_snrnaseq/run_analysis.sbatch"
