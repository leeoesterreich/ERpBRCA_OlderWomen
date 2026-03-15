#!/bin/bash
#SBATCH --job-name=spatial_biopsies
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/06_spatial_biopsies/logs/spatial_biopsies_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/06_spatial_biopsies/logs/spatial_biopsies_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Spatial Biopsies Analysis Pipeline
# Runs immune secretion spatial plots and immune pathway enrichment

set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/06_spatial_biopsies"
cd "$SCRIPT_DIR"

# Create directories
mkdir -p logs figures/immune_secretion
mkdir -p figures/immune_pathways/{enrichment,prerank,dotplots}

echo "=============================================="
echo "Spatial Biopsies Analysis Pipeline"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Working directory: $SCRIPT_DIR"
echo "Start time: $(date)"
echo ""

# Load conda and activate environment
set +u
source ~/.bashrc
eval "$(conda shell.bash hook)"
set -u
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

echo "Python version: $(python --version)"
echo ""

# Verify data symlink
if [ ! -f data/biopsy_adatas.pkl ]; then
    echo "ERROR: data/biopsy_adatas.pkl not found (symlink broken?)"
    exit 1
fi

# Run analysis scripts in order
echo "[1/2] Immune secretion spatial plots..."
python 01_immune_secretion.py 2>&1

echo "[2/2] Immune pathway enrichment..."
python 02_immune_pathways.py "$@" 2>&1

# Write completion marker
MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
touch "$MARKER_DIR/06_spatial_biopsies.complete"
echo "Marker written: $MARKER_DIR/06_spatial_biopsies.complete"

echo ""
echo "=============================================="
echo "Spatial Biopsies Analysis Complete"
echo "=============================================="
echo "End time: $(date)"
echo "Figures: $SCRIPT_DIR/figures/"
