#!/bin/bash
# analysis/06_spatial_biopsies/run_analysis.sh
# Runs spatial biopsies analysis pipeline
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Spatial Biopsies Analysis Pipeline ==="
echo "Working directory: $SCRIPT_DIR"
echo ""

# Create directories
mkdir -p logs figures/immune_secretion
mkdir -p figures/immune_pathways/{enrichment,prerank,dotplots}

# Verify data
if [ ! -f data/biopsy_adatas.pkl ]; then
    echo "ERROR: data/biopsy_adatas.pkl not found"
    exit 1
fi

echo "[1/2] Immune secretion spatial plots..."
python 01_immune_secretion.py

echo "[2/2] Immune pathway enrichment..."
python 02_immune_pathways.py "$@"

echo ""
echo "=== Spatial Biopsies Analysis Complete ==="
echo "Figures: $SCRIPT_DIR/figures/"
