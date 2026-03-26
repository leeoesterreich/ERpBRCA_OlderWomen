#!/bin/bash
#SBATCH --job-name=rat_snrnaseq
#SBATCH --output=logs/snrnaseq_%j.out
#SBATCH --error=logs/snrnaseq_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Rat snRNA-seq Analysis Pipeline
# Runs QC, integration, clustering, DE, and DA analyses

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Create directories
mkdir -p logs outputs

echo "=============================================="
echo "Rat snRNA-seq Analysis Pipeline"
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
conda activate erp_snrnaseq || {
    echo "Conda environment not found, attempting to create..."
    conda env create -f ../../envs/snrnaseq.yml
    conda activate erp_snrnaseq
}

echo "R version: $(R --version | head -1)"
echo ""

# Function to run step with error handling
run_step() {
    local step_num="$1"
    local script_name="$2"
    local description="$3"

    echo "----------------------------------------------"
    echo "Step $step_num: $description"
    echo "----------------------------------------------"

    if [[ ! -f "$SCRIPT_DIR/$script_name" ]]; then
        echo "ERROR: Script not found: $script_name"
        return 1
    fi

    Rscript "$SCRIPT_DIR/$script_name" 2>&1

    echo "Step $step_num completed."
    echo ""
}

# Run all analysis steps
run_step "01" "01_qc_filter.R" "QC, filtering, and doublet detection"
run_step "02" "02_normalize_integrate.R" "SCTransform normalization and Harmony integration"
run_step "03" "03_cluster_annotate.R" "UMAP, clustering, and cell type annotation"
run_step "08" "08_differential_expression.R" "Differential expression: Young vs Aged"
run_step "09" "09_differential_abundance.R" "Differential abundance testing"

echo "=============================================="
echo "Pipeline Complete"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output locations:"
echo "  - Intermediate files: $SCRIPT_DIR/outputs/"
echo "  - Results: ../../results/corrected/rat_snrnaseq/"
echo "  - Figures: ../../figures/by_analysis/rat_snrnaseq/"
