#!/bin/bash
# 00_setup_env.sh — Install R packages into erp_brca_aging conda env
# Adapted from Jian_MICA/scripts/00_setup_env.sh for monorepo
set -euo pipefail

ENV_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging"
MICA_PKG="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/MICA"

conda run -p "${ENV_PATH}" R --no-save <<'RSCRIPT'
# Install CRAN packages
cran_pkgs <- c("tidyverse", "devtools", "openxlsx", "ggbeeswarm")
for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        install.packages(pkg, repos = "https://cloud.r-project.org")
    } else {
        cat(pkg, "already installed.\n")
    }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_pkgs <- c("edgeR")
for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
        cat(pkg, "already installed.\n")
    }
}

cat("\n--- Verification ---\n")
all_needed <- c("tidyverse", "GSVA", "edgeR", "ComplexHeatmap", "circlize",
                "ggplot2", "patchwork", "ggbeeswarm", "openxlsx",
                "data.table", "devtools", "dplyr", "tidyr", "jsonlite")
for (pkg in all_needed) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    cat(sprintf("  %-20s %s\n", pkg, ifelse(ok, "OK", "MISSING")))
}
RSCRIPT

# Install MICA from local source
conda run -p "${ENV_PATH}" R --no-save -e "devtools::install('${MICA_PKG}', upgrade = 'never')"

# Verify MICA loads
conda run -p "${ENV_PATH}" R --no-save -e "library(MICA); cat('MICA loaded successfully.\n')"

echo "Environment setup complete."
