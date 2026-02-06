#!/bin/bash
#SBATCH --job-name=setup_envs
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/logs/setup_envs_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/logs/setup_envs_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Setup conda environments for ERp BRCA analyses

set -eo pipefail

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen"
cd "$SCRIPT_DIR"

mkdir -p "$SCRIPT_DIR/logs"

echo "=============================================="
echo "Setting up conda environments"
echo "=============================================="
echo "Start time: $(date)"
echo ""

# Initialize conda - handle potential unbound variables in profile scripts
set +u
source ~/.bashrc
eval "$(conda shell.bash hook)"
set -u

# Prefer mamba for faster solving
CONDA_CMD="mamba"
if ! command -v mamba &>/dev/null; then
    CONDA_CMD="conda"
fi
echo "Using: $CONDA_CMD"

# Create bulk RNA-seq environment
if conda env list | grep -q "erp_bulk_rnaseq"; then
    echo "erp_bulk_rnaseq environment already exists"
else
    echo "Creating erp_bulk_rnaseq environment..."
    $CONDA_CMD env create -f envs/bulk_rnaseq.yml
    echo "erp_bulk_rnaseq created successfully"

    # Install additional R packages
    echo "Installing additional R packages for bulk RNA-seq..."
    conda activate erp_bulk_rnaseq
    Rscript -e '
        options(repos = c(CRAN = "https://cloud.r-project.org"))
        if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(c("progeny", "qusage"), update = FALSE, ask = FALSE)
        install.packages("msigdbr", quiet = TRUE)
    '
    conda deactivate
fi

echo ""

# Create snRNA-seq environment
if conda env list | grep -q "erp_snrnaseq"; then
    echo "erp_snrnaseq environment already exists"
else
    echo "Creating erp_snrnaseq environment..."
    $CONDA_CMD env create -f envs/snrnaseq.yml
    echo "erp_snrnaseq created successfully"

    # Install additional R packages
    echo "Installing additional R packages for snRNA-seq..."
    conda activate erp_snrnaseq
    Rscript -e '
        options(repos = c(CRAN = "https://cloud.r-project.org"))
        if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install("SingleR", update = FALSE, ask = FALSE)
        install.packages(c("dittoSeq", "speckle"), quiet = TRUE)
        if (!require("remotes", quietly = TRUE)) install.packages("remotes")
        remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", quiet = TRUE)
    '
    conda deactivate
fi

echo ""
echo "=============================================="
echo "Environment setup complete"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "To run the analyses, submit:"
echo "  sbatch analysis/01_human_bulk_rnaseq/run_analysis.sbatch"
echo "  sbatch analysis/02_rat_snrnaseq/run_analysis.sbatch"
