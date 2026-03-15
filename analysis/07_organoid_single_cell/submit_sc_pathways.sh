#!/bin/bash
#SBATCH --job-name=sc_pathways
#SBATCH --output=logs/sc_pathways_%j.out
#SBATCH --error=logs/sc_pathways_%j.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set +u
source ~/.bashrc
eval "$(conda shell.bash hook)"
set -u
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/07_organoid_single_cell

python 07_single_cell_pathways.py
