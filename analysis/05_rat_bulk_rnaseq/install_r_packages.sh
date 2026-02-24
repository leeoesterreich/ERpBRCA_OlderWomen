#!/bin/bash
#SBATCH --job-name=install_r_pkgs
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/install_r_%j.out
#SBATCH --error=logs/install_r_%j.err

# Install R packages for rat bulk RNA-seq pipeline
set -eo pipefail

ENV_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging"

echo "=== Installing R packages into erp_brca_aging ==="
echo "Environment: $ENV_PATH"

# Activate conda (initialize first)
eval "$(/ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/conda shell.bash hook)"
conda activate "$ENV_PATH"

echo ""
echo "Installing via conda..."

# Install Bioconductor packages via conda (faster and more reliable)
mamba install -y -c conda-forge -c bioconda \
    bioconductor-deseq2 \
    bioconductor-biomart \
    bioconductor-genefu \
    r-pheatmap

echo ""
echo "=== Verifying installations ==="

Rscript -e "
suppressPackageStartupMessages({
  library(DESeq2)
  library(biomaRt)
  library(genefu)
  library(pheatmap)
  library(data.table)
  library(dplyr)
})
cat('All packages loaded successfully!\n')
cat('DESeq2 version:', as.character(packageVersion('DESeq2')), '\n')
cat('biomaRt version:', as.character(packageVersion('biomaRt')), '\n')
cat('genefu version:', as.character(packageVersion('genefu')), '\n')
"

echo ""
echo "=== Installation complete ==="
