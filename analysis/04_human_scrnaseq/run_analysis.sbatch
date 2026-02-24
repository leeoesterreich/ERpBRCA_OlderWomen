#!/bin/bash
#SBATCH --job-name=human_scrnaseq
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/analysis_%j.out
#SBATCH --error=logs/analysis_%j.err

# Human scRNA-seq analysis pipeline
# Requires: Raw data downloaded via 00a_download_geo.sh

# Load conda environment (before strict mode to avoid bashrc variable issues)
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

# Now enable strict mode
set -eo pipefail

# Use SLURM_SUBMIT_DIR to get the directory where sbatch was called
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
echo "Working directory: $SCRIPT_DIR"

mkdir -p logs outputs outputs/preprocessing

# Check if preprocessed data already exists
PREPROCESS_FILE="../../data/human_scrnaseq/SeuratObj_GSE176078_ERpos_AfterQCSCT.rds"
if [ -f "$PREPROCESS_FILE" ]; then
    echo "=== Step 0b: Preprocessed data exists, skipping ==="
else
    echo "=== Step 0b: Preprocess raw data ==="
    Rscript 00b_preprocess_seurat.R
fi

echo "=== Step 1: Load and subset data ==="
Rscript 01_load_subset_data.R

# Note: Harmony integration removed (not used per methods)

echo "=== Step 3: Cell type annotation ==="
Rscript 03_cell_type_annotation.R

echo "=== Step 4: Cell fractions ==="
Rscript 04_cell_fractions.R

echo "=== Step 5: Gene expression ==="
Rscript 05_gene_expression_violin.R

echo "=== Step 6: GSVA ==="
Rscript 06_run_gsva.R

echo "=== Step 7: PROGENy ==="
Rscript 07_run_progeny.R

echo "=== Step 8: WCSEA ==="
Rscript 08_run_wcsea.R

echo "=== Step 9: CellPhoneDB prep ==="
Rscript 09_cellphonedb_prep.R

echo "=== Analysis complete ==="
