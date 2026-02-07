#!/bin/bash
#SBATCH --job-name=comparison_report
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/03_comparison/logs/comparison_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/03_comparison/logs/comparison_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Generate comparison report: Original vs Corrected results

set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/03_comparison"
cd "$SCRIPT_DIR"

mkdir -p logs

echo "=============================================="
echo "Generating Comparison Report"
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

Rscript "$SCRIPT_DIR/01_generate_comparison_report.R" 2>&1

echo ""
echo "=============================================="
echo "Complete"
echo "=============================================="
echo "End time: $(date)"
