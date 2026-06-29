#!/bin/bash
#SBATCH --job-name=scrnaseq_preproc
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/preprocessing_%j.out
#SBATCH --error=logs/preprocessing_%j.err

# Stage 1: Load, annotate, cell fractions, gene expression
# Produces: outputs/seurat_annotated.rds (checkpoint for Stage 2)

# Load conda environment on HTC cluster
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin:$PATH"
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging

# Now enable strict mode
set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/04_human_scrnaseq"
cd "$SCRIPT_DIR"
echo "Working directory: $SCRIPT_DIR"

mkdir -p logs outputs outputs/preprocessing figures

# Step 0b: Preprocess raw data (skip if already done)
PREPROCESS_FILE="../../data/human_scrnaseq/SeuratObj_Xu2024_HRpos_AfterQCSCT.rds"
if [ -f "$PREPROCESS_FILE" ]; then
    echo "=== Step 0b: Preprocessed data exists, skipping ==="
else
    echo "=== Step 0b: Preprocess raw data ==="
    Rscript 00b_preprocess_seurat.R
fi

echo "=== Step 1: Load and subset data ==="
Rscript 01_load_subset_data.R

# Note: Harmony integration removed (not used per methods)

echo "=== Step 3: Cell type annotation ==="
Rscript 03_cell_type_annotation.R

echo "=== Step 4: Cell fractions ==="
Rscript 04_cell_fractions.R

echo "=== Step 5: Gene expression ==="
Rscript 05_gene_expression_violin.R

# Validate checkpoint file was created
if [ ! -f "outputs/seurat_annotated.rds" ]; then
    echo "ERROR: seurat_annotated.rds not created — Stage 1 failed"
    exit 1
fi

# Write preprocessing completion marker
MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
echo "Stage 1 (preprocessing) Job ${SLURM_JOB_ID:-local} completed at $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$MARKER_DIR/04_human_scrnaseq_preprocessing.complete"
echo "Marker written: $MARKER_DIR/04_human_scrnaseq_preprocessing.complete"

echo "=== Preprocessing complete ==="
