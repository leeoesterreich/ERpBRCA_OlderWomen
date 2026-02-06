#!/bin/bash
# analysis/01_human_bulk_rnaseq/run_analysis.sh
# Runs complete human bulk RNA-seq analysis pipeline
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Human Bulk RNA-seq Analysis Pipeline ==="
echo "Working directory: $SCRIPT_DIR"
echo ""

# Create outputs directory
mkdir -p outputs

# Run analysis scripts in order
echo "[1/5] Preprocessing..."
Rscript 01_preprocess.R

echo "[2/5] Running GSVA..."
Rscript 02_run_gsva.R

echo "[3/5] Running PROGENy..."
Rscript 03_run_progeny.R

echo "[4/5] Calculating correlations..."
Rscript 04_correlations.R

echo "[5/5] Generating visualizations..."
Rscript 05_visualize.R

echo ""
echo "=== Human Bulk RNA-seq Analysis Complete ==="
echo "Outputs: $SCRIPT_DIR/outputs/"
echo "Figures: ../../figures/by_analysis/human_bulk_rnaseq/"
