#!/bin/bash
#SBATCH --job-name=cpdb_v4
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/cpdb_v4_%j.out
#SBATCH --error=logs/cpdb_v4_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

source ~/.bashrc
conda activate cellphonedb_v4

set -eo pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname $(readlink -f $0))}"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "CellPhoneDB v4 Analysis"
echo "=========================================="
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "Python: $(which python)"
python -c "import pkg_resources; print(f'CellPhoneDB version: {pkg_resources.get_distribution(\"cellphonedb\").version}')"

# Run CellPhoneDB v4
python 16b_run_cellphonedb_v4.py

echo ""
echo "Listing output files:"
ls -la outputs/cellphonedb_v4/results_*/ 2>/dev/null || echo "No output files yet"

echo ""
echo "Done: $(date)"
