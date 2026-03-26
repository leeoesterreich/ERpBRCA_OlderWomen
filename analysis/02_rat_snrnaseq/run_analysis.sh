#!/bin/bash
# analysis/02_rat_snrnaseq/run_analysis.sh
# Master script to run the complete rat snRNA-seq analysis pipeline
#
# This script runs all analysis steps in sequence:
# 1. QC and filtering with doublet detection
# 2. Normalization and Harmony integration
# 3. Clustering and cell type annotation
# 4. (Steps 4-7 reserved for future analyses)
# 8. Differential expression (Young vs Aged)
# 9. Differential abundance testing
#
# Usage:
#   ./run_analysis.sh           # Run all steps
#   ./run_analysis.sh 01        # Run only step 01
#   ./run_analysis.sh 01 03     # Run steps 01 through 03

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "=============================================="
echo "Rat snRNA-seq Analysis Pipeline"
echo "=============================================="
echo "Script directory: $SCRIPT_DIR"
echo "Project root: $PROJECT_ROOT"
echo "Start time: $(date)"
echo ""

# Function to run an R script with error handling
run_step() {
    local step_num="$1"
    local script_name="$2"
    local description="$3"

    echo "----------------------------------------------"
    echo "Step $step_num: $description"
    echo "----------------------------------------------"

    local script_path="$SCRIPT_DIR/$script_name"

    if [[ ! -f "$script_path" ]]; then
        echo "ERROR: Script not found: $script_path"
        return 1
    fi

    echo "Running: $script_path"
    echo ""

    if Rscript "$script_path"; then
        echo ""
        echo "Step $step_num completed successfully."
    else
        echo ""
        echo "ERROR: Step $step_num failed!"
        return 1
    fi

    echo ""
}

# Parse arguments for step selection
START_STEP="${1:-01}"
END_STEP="${2:-09}"

echo "Running steps $START_STEP through $END_STEP"
echo ""

# Step 1: QC and filtering with doublet detection
if [[ "$START_STEP" -le 1 && "$END_STEP" -ge 1 ]]; then
    run_step "01" "01_qc_filter.R" "QC, filtering, and doublet detection"
fi

# Step 2: Normalization and integration
if [[ "$START_STEP" -le 2 && "$END_STEP" -ge 2 ]]; then
    run_step "02" "02_normalize_integrate.R" "SCTransform normalization and Harmony integration"
fi

# Step 3: Clustering and annotation
if [[ "$START_STEP" -le 3 && "$END_STEP" -ge 3 ]]; then
    run_step "03" "03_cluster_annotate.R" "UMAP, clustering, and cell type annotation"
fi

# Steps 4-7 are reserved for future analyses
# (e.g., trajectory analysis, gene regulatory networks, etc.)

# Step 8: Differential expression
if [[ "$START_STEP" -le 8 && "$END_STEP" -ge 8 ]]; then
    run_step "08" "08_differential_expression.R" "Differential expression: Young vs Aged"
fi

# Step 9: Differential abundance
if [[ "$START_STEP" -le 9 && "$END_STEP" -ge 9 ]]; then
    run_step "09" "09_differential_abundance.R" "Differential abundance testing"
fi

echo "=============================================="
echo "Pipeline Complete"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output locations:"
echo "  - Intermediate files: $SCRIPT_DIR/outputs/"
echo "  - Results: $PROJECT_ROOT/results/corrected/rat_snrnaseq/"
echo "  - Figures: $PROJECT_ROOT/figures/by_analysis/rat_snrnaseq/"
echo ""
echo "Key output files:"
echo "  - seurat_qc_filtered.rds (post-QC, doublets removed)"
echo "  - seurat_integrated.rds (normalized, integrated)"
echo "  - seurat_annotated.rds (clustered, annotated)"
echo "  - DE_results_by_celltype.csv (differential expression)"
echo "  - DA_results_celltypes.csv (differential abundance)"
