#!/bin/bash
# Master script for manuscript validation pipeline

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

cd "$PROJ_DIR"

echo "========================================"
echo "Manuscript Validation Pipeline"
echo "========================================"
echo "Project: $PROJ_DIR"
echo "Started: $(date)"
echo ""

echo "=== Step 1: Extract figures from PowerPoint ==="
python3 scripts/validation/01_extract_figures.py
echo ""

echo "=== Step 2: Extract legends and claims from DOCX ==="
python3 scripts/validation/02_extract_claims.py
echo ""

echo "=== Step 3: Check output status ==="
python3 scripts/validation/02b_check_outputs.py
echo ""

echo "=== Step 4: Numerical comparison (Tier 1) ==="
python3 scripts/validation/03_numerical_comparison.py
echo ""

echo "=== Step 5: Visual comparison (Tier 2) ==="
python3 scripts/validation/04_visual_comparison.py
echo ""

echo "=== Step 6: Claim verification (Tier 3) ==="
python3 scripts/validation/05_claim_verification.py
echo ""

echo "=== Step 7: Generate report ==="
python3 scripts/validation/06_generate_report.py
echo ""

echo "========================================"
echo "Validation Complete"
echo "========================================"
echo "Report: $PROJ_DIR/validation/report/validation_report.md"
echo "Finished: $(date)"
