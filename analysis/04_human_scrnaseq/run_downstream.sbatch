#!/bin/bash
#SBATCH --job-name=scrnaseq_downstream
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/downstream_%j.out
#SBATCH --error=logs/downstream_%j.err

# Stage 2: GSVA, PROGENy, WCSEA, CellPhoneDB prep, TAM export
# Depends on: outputs/seurat_annotated.rds from Stage 1

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

# Validate Stage 1 checkpoint exists
if [ ! -f "outputs/seurat_annotated.rds" ]; then
    echo "ERROR: outputs/seurat_annotated.rds not found. Run Stage 1 (run_preprocessing.sbatch) first."
    exit 1
fi

echo "=== Step 6: GSVA (pseudobulk only) ==="
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

  if (!file.exists("outputs/seurat_annotated.rds")) {
    stop("seurat_annotated.rds not found")
  }

  seurat <- readRDS("outputs/seurat_annotated.rds")

  if (!"CellTypeAnnot" %in% colnames(seurat@meta.data)) {
    stop("Column CellTypeAnnot not found in Seurat metadata")
  }

  mac_cells <- WhichCells(seurat, expression = CellTypeAnnot == "Macrophage")

  if (length(mac_cells) == 0) {
    warning("No macrophage cells found, creating empty TAM analysis")
    tam_metrics <- data.frame(sample = character(), age_group = character(), n_cells = integer())
  } else {
    mac_subset <- subset(seurat, cells = mac_cells)

    tam_metrics <- mac_subset@meta.data %>%
      group_by(orig.ident) %>%
      summarise(
        n_cells = n(),
        .groups = "drop"
      )

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

MARKER_DIR="$(dirname "$SCRIPT_DIR")/.pipeline_markers"
mkdir -p "$MARKER_DIR"
echo "Stage 2 (downstream) Job ${SLURM_JOB_ID:-local} completed at $(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$MARKER_DIR/04_human_scrnaseq.complete"
echo "Marker written: $MARKER_DIR/04_human_scrnaseq.complete"

echo "=== Downstream analysis complete ==="
