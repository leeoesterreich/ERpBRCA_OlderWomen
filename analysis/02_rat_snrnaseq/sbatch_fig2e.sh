#!/usr/bin/env bash
#SBATCH --job-name=rat_fig2e
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/rat_fig2e_%j.out
#SBATCH --error=logs/rat_fig2e_%j.err

set -eo pipefail
unset LD_LIBRARY_PATH
source ${CONDA_PREFIX_ROOT:-$HOME/miniconda3}/etc/profile.d/conda.sh
conda activate ${ERP_BRCA_ENV:-erp_brca_aging}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/../.."
mkdir -p logs

echo "Working dir: $(pwd)"
echo "R: $(which Rscript)"
echo "Start: $(date)"

Rscript analysis/02_rat_snrnaseq/10_multicelltype_pathway.R

echo "End: $(date)"
