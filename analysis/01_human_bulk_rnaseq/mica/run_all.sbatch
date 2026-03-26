#!/bin/bash
#SBATCH --job-name=mica_analysis
#SBATCH --partition=htc
#SBATCH --time=4:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/mica_%j.out
#SBATCH --error=logs/mica_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

# Activate environment
source ~/.bashrc
conda activate erp_brca_aging

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
cd "$SCRIPT_DIR"
mkdir -p logs outputs figures

echo "=== MICA Analysis Pipeline ==="
echo "Start time: $(date)"
echo "Working directory: $(pwd)"

# Step 1: Prepare data and compute GSVA
echo "Step 1: Preparing data and computing GSVA..."
Rscript 01_prep_data.R

# Step 2: GSVA (integrated into step 1)
echo "Step 2: GSVA (integrated into 01_prep_data.R)..."
Rscript 02_run_gsva.R

# Step 3: Run MICA concordance analysis
echo "Step 3: Running MICA..."
Rscript 03_run_mica.R

# Step 4: Generate Figure 4 (estrogen/HSD17B)
echo "Step 4: Generating Figure 4..."
Rscript 04_fig4_estrogen.R

# Step 5: Generate Figure 5 (immune/inflammation)
echo "Step 5: Generating Figure 5..."
Rscript 05_fig5_immune.R

# Step 6: Extract manuscript claims
echo "Step 6: Extracting manuscript claims..."
python 06_extract_claims.py

# Step 7: Statistical audit
echo "Step 7: Running statistical audit..."
Rscript 07_statistical_audit.R

# Step 8: Compile final report
echo "Step 8: Compiling report..."
python 08_compile_report.py

# Write completion marker
MARKER_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)/.pipeline_markers"
mkdir -p "$MARKER_DIR"
date > "$MARKER_DIR/01_mica.complete"
echo "Completion marker written to $MARKER_DIR/01_mica.complete"

echo "=== Pipeline complete ==="
echo "End time: $(date)"
