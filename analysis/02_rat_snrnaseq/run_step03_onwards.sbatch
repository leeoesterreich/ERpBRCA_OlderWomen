#!/bin/bash
#SBATCH --job-name=rat_snrnaseq_step03
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/02_rat_snrnaseq/logs/snrnaseq_step03_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/02_rat_snrnaseq/logs/snrnaseq_step03_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Rat snRNA-seq Analysis Pipeline - Step 03 onwards
# Skips QC and normalization since they're already complete

set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/02_rat_snrnaseq"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "Rat snRNA-seq Analysis Pipeline (Step 03+)"
echo "=============================================="
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start time: $(date)"
echo ""

# Load conda and activate environment
set +u
source ~/.bashrc
eval "$(conda shell.bash hook)"
set -u
conda activate erp_snrnaseq

echo "R version: $(R --version | head -1)"
echo ""

# Install missing packages
echo "Installing missing R packages..."
Rscript -e "
if (!require('pheatmap', quietly = TRUE)) install.packages('pheatmap', repos='https://cloud.r-project.org')
if (!require('ggrepel', quietly = TRUE)) install.packages('ggrepel', repos='https://cloud.r-project.org')
if (!require('speckle', quietly = TRUE)) BiocManager::install('speckle', ask=FALSE, update=FALSE)
cat('Package installation complete\n')
" 2>&1

# Verify step 02 output exists
if [[ ! -f "$SCRIPT_DIR/outputs/seurat_integrated.rds" ]]; then
    echo "ERROR: Step 02 output not found. Run steps 01-02 first."
    exit 1
fi

run_step() {
    local step_num="$1"
    local script_name="$2"
    local description="$3"

    echo "----------------------------------------------"
    echo "Step $step_num: $description"
    echo "----------------------------------------------"

    Rscript "$SCRIPT_DIR/$script_name" 2>&1

    echo "Step $step_num completed."
    echo ""
}

# Run steps 03, 08, 09
run_step "03" "03_cluster_annotate.R" "UMAP, clustering, and cell type annotation"
run_step "08" "08_differential_expression.R" "Differential expression: Young vs Aged"
run_step "09" "09_differential_abundance.R" "Differential abundance testing"

echo "=============================================="
echo "Pipeline Complete"
echo "=============================================="
echo "End time: $(date)"
