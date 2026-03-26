#!/bin/bash
# scripts/run_all.sh
# Master script to run complete analysis pipeline
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "=============================================="
echo "  ERpBRCA OlderWomen Analysis Pipeline"
echo "=============================================="
echo "Project root: $PROJECT_ROOT"
echo "Started at: $(date)"
echo ""

# Parse arguments
RUN_DOWNLOAD=true
RUN_BULK=true
RUN_SNRNA=true
RUN_WES=true
RUN_ORGANOID=true
RUN_COMPARE=true

while [[ $# -gt 0 ]]; do
  case $1 in
    --skip-download) RUN_DOWNLOAD=false; shift ;;
    --bulk-only) RUN_SNRNA=false; RUN_WES=false; RUN_ORGANOID=false; shift ;;
    --snrna-only) RUN_BULK=false; RUN_WES=false; RUN_ORGANOID=false; shift ;;
    --skip-organoid) RUN_ORGANOID=false; shift ;;
    --skip-compare) RUN_COMPARE=false; shift ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# Step 1: Download data
if [ "$RUN_DOWNLOAD" = true ]; then
  echo ""
  echo "[1/7] Downloading data from GEO..."
  echo "----------------------------------------------"
  bash data/download_all.sh
else
  echo "[1/7] Skipping data download..."
fi

# Step 2: Human bulk RNA-seq
if [ "$RUN_BULK" = true ]; then
  echo ""
  echo "[2/7] Running human bulk RNA-seq analysis..."
  echo "----------------------------------------------"
  cd analysis/01_human_bulk_rnaseq
  bash run_analysis.sh
  cd "$PROJECT_ROOT"
else
  echo "[2/7] Skipping bulk RNA-seq analysis..."
fi

# Step 3: Rat snRNA-seq
if [ "$RUN_SNRNA" = true ]; then
  echo ""
  echo "[3/7] Running rat snRNA-seq analysis..."
  echo "----------------------------------------------"
  cd analysis/02_rat_snrnaseq
  bash run_analysis.sh
  cd "$PROJECT_ROOT"
else
  echo "[3/7] Skipping snRNA-seq analysis..."
fi

# Step 4: Rat WES
if [ "$RUN_WES" = true ]; then
  echo ""
  echo "[4/7] Running rat WES analysis..."
  echo "----------------------------------------------"
  cd analysis/03_rat_wes
  bash run_analysis.sh
  cd "$PROJECT_ROOT"
else
  echo "[4/7] Skipping WES analysis..."
fi

# Step 5: Rat bulk RNA-seq (placeholder — uses same structure as other sections)
echo ""
echo "[5/7] Rat bulk RNA-seq (run separately if needed)..."
echo "----------------------------------------------"

# Step 6: Organoid single-cell
if [ "$RUN_ORGANOID" = true ]; then
  echo ""
  echo "[6/7] Running organoid single-cell analysis..."
  echo "----------------------------------------------"
  cd analysis/07_organoid_single_cell
  bash run_analysis.sh
  cd "$PROJECT_ROOT"
else
  echo "[6/7] Skipping organoid single-cell analysis..."
fi

# Step 7: Generate comparison report
if [ "$RUN_COMPARE" = true ]; then
  echo ""
  echo "[7/7] Generating comparison report..."
  echo "----------------------------------------------"
  Rscript scripts/compare_results.R
else
  echo "[7/7] Skipping comparison report..."
fi

# Organize figures
echo ""
echo "Organizing figures for manuscript..."
echo "----------------------------------------------"
Rscript scripts/organize_figures.R 2>/dev/null || echo "Figure organization script not yet implemented"

echo ""
echo "=============================================="
echo "  Pipeline Complete"
echo "=============================================="
echo "Finished at: $(date)"
echo ""
echo "Results:"
echo "  - Comparison report: results/comparison/comparison_report.pdf"
echo "  - Figures: figures/manuscript/"
echo "  - Corrected results: results/corrected/"
echo ""
