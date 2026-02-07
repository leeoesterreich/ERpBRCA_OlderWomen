#!/bin/bash
#SBATCH --job-name=rat_step08
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/02_rat_snrnaseq/logs/snrnaseq_step08_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/02_rat_snrnaseq/logs/snrnaseq_step08_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Re-run step 08 with future.globals.maxSize fix for Luminal cells

set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/02_rat_snrnaseq"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "Rat snRNA-seq Step 08 Re-run"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo ""

# Load conda and activate environment
set +u
source ~/.bashrc
eval "$(conda shell.bash hook)"
set -u
conda activate erp_snrnaseq

echo "R version: $(R --version | head -1)"
echo ""

echo "----------------------------------------------"
echo "Step 08: Differential expression with maxSize fix"
echo "----------------------------------------------"

Rscript "$SCRIPT_DIR/08_differential_expression.R" 2>&1

echo "Step 08 completed."
echo ""

echo "=============================================="
echo "Complete"
echo "=============================================="
echo "End time: $(date)"
