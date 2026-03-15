#!/bin/bash
#SBATCH --job-name=neil_scrna_analysis
#SBATCH --output=logs/analysis_%j.out
#SBATCH --error=logs/analysis_%j.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Activate environment
set +u
source ~/.bashrc
eval "$(conda shell.bash hook)"
set -u
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/07_organoid_single_cell

echo "Starting analysis pipeline at $(date)"
# echo "Step 1: QC"
# python 02_qc.py

# echo "Step 2: Preprocessing"
# python 03_preprocess.py

echo "Step 3: Single-cell pathway analysis"
python 07_single_cell_pathways.py

echo "Pipeline complete at $(date)!"
