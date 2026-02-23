#!/bin/bash
#SBATCH --job-name=rat_bulk_R
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/R_analysis_%j.out
#SBATCH --error=logs/R_analysis_%j.err

# Rat bulk RNA-seq R analysis
set -euo pipefail

module purge
module load r/4.2.0

# Use absolute path for SLURM compatibility
SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/05_rat_bulk_rnaseq"
cd "$SCRIPT_DIR"

mkdir -p logs

echo "=== DESeq2 Analysis ==="
Rscript 03_deseq2.R

echo "=== PAM50 Subtyping ==="
Rscript 04_pam50_subtyping.R

echo "=== Validation ==="
Rscript 05_validate_vs_rahul.R

echo "=== Analysis complete ==="
