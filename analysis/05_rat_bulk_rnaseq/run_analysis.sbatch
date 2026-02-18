#!/bin/bash
#SBATCH --job-name=rat_bulk_analysis
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/analysis_%j.out
#SBATCH --error=logs/analysis_%j.err

# Rat bulk RNA-seq R analysis
set -euo pipefail

module purge
module load r/4.2.0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p logs

echo "=== DESeq2 Analysis ==="
Rscript 03_deseq2.R

echo "=== PAM50 Subtyping ==="
Rscript 04_pam50_subtyping.R

echo "=== Analysis complete ==="
