#!/bin/bash
#SBATCH --job-name=human_bulk_rnaseq
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/01_human_bulk_rnaseq/logs/bulk_rnaseq_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/01_human_bulk_rnaseq/logs/bulk_rnaseq_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Human Bulk RNA-seq Analysis Pipeline
# Runs preprocessing, GSVA, PROGENy, correlations, and visualization

set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/01_human_bulk_rnaseq"
cd "$SCRIPT_DIR"

# Create log directory
mkdir -p logs outputs

echo "=============================================="
echo "Human Bulk RNA-seq Analysis Pipeline"
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
conda activate erp_bulk_rnaseq || {
    echo "Conda environment not found, attempting to create..."
    conda env create -f ../../envs/bulk_rnaseq.yml
    conda activate erp_bulk_rnaseq
}

echo "R version: $(R --version | head -1)"
echo ""

# Run analysis scripts in order
echo "[1/5] Preprocessing..."
Rscript 01_preprocess.R 2>&1

echo "[2/5] Running GSVA..."
Rscript 02_run_gsva.R 2>&1

echo "[3/5] Running PROGENy..."
Rscript 03_run_progeny.R 2>&1

echo "[4/5] Calculating correlations..."
Rscript 04_correlations.R 2>&1

echo "[5/5] Generating visualizations..."
Rscript 05_visualize.R 2>&1

echo "[6/6] Exporting CSVs for validation..."
Rscript -e '
# Export correlation results as CSV
cor_data <- readRDS("outputs/correlation_results.rds")
write.csv(cor_data, "outputs/correlation_results.csv", row.names = FALSE)

# Export GSVA results as CSV (for MICA analysis)
gsva_data <- readRDS("outputs/gsva_estrogen_pathways.rds")
gsva_df <- as.data.frame(t(gsva_data))
gsva_df$Sample <- rownames(gsva_df)
write.csv(gsva_df, "outputs/mica_results.csv", row.names = FALSE)

# Export HSD17B7 correlations (subset of correlation_results)
hsd_cor <- cor_data[cor_data$GeneSymb == "HSD17B7", ]
write.csv(hsd_cor, "outputs/hsd17b7_correlations.csv", row.names = FALSE)

cat("CSV exports complete\n")
'

# Write completion marker
MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
touch "$MARKER_DIR/01_human_bulk_rnaseq.complete"
echo "Marker written: $MARKER_DIR/01_human_bulk_rnaseq.complete"

echo ""
echo "=============================================="
echo "Human Bulk RNA-seq Analysis Complete"
echo "=============================================="
echo "End time: $(date)"
echo "Outputs: $SCRIPT_DIR/outputs/"
