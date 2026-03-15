#!/bin/bash
# analysis/07_organoid_single_cell/run_analysis.sh
# Runs organoid single-cell analysis pipeline (scripts 02-11)
# Cell Ranger (01) and reference building must be run separately first.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Organoid Single-Cell Analysis Pipeline ==="
echo "Working directory: $SCRIPT_DIR"
echo ""

# Create directories
mkdir -p outputs/{qc,de,pathways,cell_cycle,single_cell_pathways,inhibitor_mechanism,heterogeneity,proliferation}
mkdir -p figures/{qc,umap,de,pathways,cell_cycle,single_cell_pathways,inhibitor_mechanism,heterogeneity,proliferation}
mkdir -p logs

echo "[1/10] QC and filtering..."
python 02_qc.py 2>&1

echo "[2/10] Preprocessing (normalize, HVG, PCA, UMAP)..."
python 03_preprocess.py 2>&1

echo "[3/10] Pseudobulk differential expression..."
python 04_pseudobulk_de.py 2>&1

echo "[4/10] GSEA and decoupler pathway analysis..."
python 05_pathways.py 2>&1

echo "[5/10] Cell cycle distribution..."
python 06_cell_cycle.py 2>&1

echo "[6/10] Single-cell pathway scoring..."
python 07_single_cell_pathways.py 2>&1

echo "[7/10] Inhibitor mechanism exploration..."
python 08_inhibitor_mechanism_exploration.py 2>&1

echo "[8/10] Pathway enrichment analysis..."
python 09_pathway_enrichment_analysis.py 2>&1

echo "[9/10] Heterogeneity analysis..."
python 10_heterogeneity_analysis.py 2>&1

echo "[10/10] Proliferation analysis..."
python 11_proliferation_analysis.py 2>&1

echo ""
echo "=== Organoid Single-Cell Analysis Complete ==="
echo "Outputs: $SCRIPT_DIR/outputs/"
echo "Figures: $SCRIPT_DIR/figures/"
