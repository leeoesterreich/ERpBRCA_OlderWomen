#!/bin/bash
#SBATCH --job-name=human_scrnaseq
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/analysis_%j.out
#SBATCH --error=logs/analysis_%j.err

# Human scRNA-seq analysis pipeline
# Requires: Raw data downloaded via 00a_download_geo.sh

# Load conda environment (before strict mode to avoid bashrc variable issues)
source ~/.bashrc
conda activate erp_brca_aging

# Now enable strict mode
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
echo "Working directory: $SCRIPT_DIR"

mkdir -p logs outputs outputs/preprocessing figures

# Check if preprocessed data already exists
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

echo "=== Step 6: GSVA (pseudobulk only — single-cell mode exceeds 256GB) ==="
Rscript 06_run_gsva.R --mode=pseudobulk

echo "=== Step 7: PROGENy ==="
Rscript 07_run_progeny.R

echo "=== Step 8: WCSEA ==="
Rscript 08_run_wcsea.R

echo "=== Step 9: CellPhoneDB prep ==="
Rscript 09_cellphonedb_prep.R

echo "=== Step 10: Export TAM analysis CSV ==="
Rscript -e '
tryCatch({
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(dplyr))

  # Check input file exists
  if (!file.exists("outputs/seurat_annotated.rds")) {
    stop("seurat_annotated.rds not found")
  }

  # Load annotated object
  seurat <- readRDS("outputs/seurat_annotated.rds")

  # Validate required column
  if (!"CellTypeAnnot" %in% colnames(seurat@meta.data)) {
    stop("Column CellTypeAnnot not found in Seurat metadata")
  }

  # Get macrophage cells
  mac_cells <- WhichCells(seurat, expression = CellTypeAnnot == "Macrophage")

  if (length(mac_cells) == 0) {
    warning("No macrophage cells found, creating empty TAM analysis")
    tam_metrics <- data.frame(sample = character(), age_group = character(), n_cells = integer())
  } else {
    mac_subset <- subset(seurat, cells = mac_cells)

    # Calculate TAM metrics per sample
    tam_metrics <- mac_subset@meta.data %>%
      group_by(orig.ident) %>%
      summarise(
        n_cells = n(),
        .groups = "drop"
      )

    # Add AgeGroup if available
    if ("AgeGroup" %in% colnames(mac_subset@meta.data)) {
      age_mapping <- mac_subset@meta.data %>%
        select(orig.ident, AgeGroup) %>%
        distinct()
      tam_metrics <- left_join(tam_metrics, age_mapping, by = "orig.ident")
    }
  }

  write.csv(tam_metrics, "outputs/tam_analysis.csv", row.names = FALSE)
  cat("TAM analysis exported:", nrow(tam_metrics), "samples\n")
}, error = function(e) {
  message(paste("ERROR in TAM export:", e$message))
  quit(status = 1)
})
' || { echo "ERROR: TAM export failed"; exit 1; }

echo "=== Step 11: Ensure PNG figures exist ==="
# Convert any PDFs that don't have PNG counterparts (fallback)
for pdf in outputs/*.pdf; do
    [ -f "$pdf" ] || continue
    base=$(basename "$pdf" .pdf)
    png="figures/${base}.png"
    if [ ! -f "$png" ]; then
        echo "Converting $pdf -> $png"
        env -i PATH=/usr/bin:/bin pdftoppm -png -r 300 -singlefile "$pdf" "figures/${base}" 2>/dev/null \
            || env -i PATH=/usr/bin:/bin convert -density 300 "$pdf" "$png" 2>/dev/null \
            || echo "WARNING: Could not convert $pdf"
    fi
done

# Write completion marker (only if we get here = success)
# Note: Single dirname to go from 04_human_scrnaseq/ up to analysis/
MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
echo "Job ${SLURM_JOB_ID:-local} completed at $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$MARKER_DIR/04_human_scrnaseq.complete"
echo "Marker written: $MARKER_DIR/04_human_scrnaseq.complete"

echo "=== Analysis complete ==="
