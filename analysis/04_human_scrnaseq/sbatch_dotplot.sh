#!/bin/bash
#SBATCH --job-name=fig7d_dotplot
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/fig7d_dotplot_%j.out
#SBATCH --error=logs/fig7d_dotplot_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname $(readlink -f $0))}"
cd "$SCRIPT_DIR"

source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

echo "Running CellPhoneDB dot plot..."
Rscript 17_cellphonedb_dotplot.R

echo "Done"
ls -la figures/fig7d*.png 2>/dev/null || echo "No output file found"
