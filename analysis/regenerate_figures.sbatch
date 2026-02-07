#!/bin/bash
#SBATCH --job-name=regen_figs
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/logs/regen_figs_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/logs/regen_figs_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Regenerate all figures in PNG and SVG formats

set -eo pipefail

PROJECT_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen"
cd "$PROJECT_ROOT"

mkdir -p analysis/logs

echo "=============================================="
echo "Regenerating Figures (PNG/SVG)"
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

# Install required packages
Rscript -e "
if (!require('patchwork', quietly = TRUE)) install.packages('patchwork', repos='https://cloud.r-project.org')
if (!require('svglite', quietly = TRUE)) install.packages('svglite', repos='https://cloud.r-project.org')
cat('Packages ready\n')
" 2>&1

echo "----------------------------------------------"
echo "Step 1: Human Bulk RNA-seq Figures"
echo "----------------------------------------------"
Rscript analysis/01_human_bulk_rnaseq/05_visualize.R 2>&1

echo ""
echo "----------------------------------------------"
echo "Step 2: Rat snRNA-seq DE Figures"
echo "----------------------------------------------"
Rscript analysis/02_rat_snrnaseq/08_differential_expression.R 2>&1

echo ""
echo "----------------------------------------------"
echo "Step 3: Rat snRNA-seq DA Figures"
echo "----------------------------------------------"
Rscript analysis/02_rat_snrnaseq/09_differential_abundance.R 2>&1

echo ""
echo "----------------------------------------------"
echo "Listing generated figures"
echo "----------------------------------------------"
echo "Human Bulk RNA-seq:"
ls -la figures/by_analysis/human_bulk_rnaseq/*.{png,svg} 2>/dev/null || echo "No PNG/SVG files yet"
echo ""
echo "Rat snRNA-seq:"
ls -la figures/by_analysis/rat_snrnaseq/*.{png,svg} 2>/dev/null || echo "No PNG/SVG files yet"

echo ""
echo "=============================================="
echo "Complete"
echo "=============================================="
echo "End time: $(date)"
