#!/usr/bin/env bash
# prepare_data_asset.sh - Package input data for CodeOcean capsule
#
# Usage:
#   ./scripts/prepare_data_asset.sh --tier 1 --output /path/to/staging
#   ./scripts/prepare_data_asset.sh --tier 2 --output /path/to/staging
#
# Tier 1 (~2.3 GB): count matrices + figure-regen CSVs
# Tier 2 (~20-35 GB): adds precomputed Seurat/h5ad intermediates
set -euo pipefail

# ── Repo root (resolve relative to this script) ────────────────────────────
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ── External data paths (set env vars or edit these defaults) ─────────────
# Xu et al. 2024 primary breast tumor atlas (barcodes.tsv, genes.tsv, matrix.mtx, metadata.csv)
XU_ATLAS="${XU_ATLAS:-/path/to/Xu_etal_Primary_Breast_Tumor_Atlas_2024}"
# Organoid scRNA-seq Cell Ranger project directory (contains pool1_flex/, pool2_flex/)
NEIL_PROJECT="${NEIL_PROJECT:-/path/to/organoid_scrnaseq_cellranger_project}"
# HCC22-088 spatial biopsy AnnData pickle
BIOPSY_PKL="${BIOPSY_PKL:-/path/to/biopsy_adatas.pkl}"

# ── CLI parsing ─────────────────────────────────────────────────────────────
TIER=""
OUTPUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tier)  TIER="$2";   shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 --tier {1|2} --output /path/to/staging"
            exit 0 ;;
        *) echo "ERROR: Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$TIER" || -z "$OUTPUT" ]]; then
    echo "ERROR: --tier and --output are required."
    echo "Usage: $0 --tier {1|2} --output /path/to/staging"
    exit 1
fi
if [[ "$TIER" != "1" && "$TIER" != "2" ]]; then
    echo "ERROR: --tier must be 1 or 2 (got: $TIER)"
    exit 1
fi

# ── Helpers ─────────────────────────────────────────────────────────────────
WARNINGS=0
SECTION_FILE_COUNTS=()

warn() { echo "  [WARN] $*"; ((WARNINGS++)) || true; }
info() { echo "  [INFO] $*"; }

# Copy a source to dest, creating parent dirs. Warn and skip if source missing.
safe_copy() {
    local src="$1" dest="$2"
    if [[ ! -e "$src" ]]; then
        warn "Source not found, skipping: $src"
        return 1
    fi
    mkdir -p "$(dirname "$dest")"
    if [[ -d "$src" ]]; then
        cp -a "$src" "$dest"
    else
        cp -a "$src" "$dest"
    fi
    return 0
}

# Copy directory contents (not the dir itself) into dest.
safe_copy_contents() {
    local src_dir="$1" dest_dir="$2"
    if [[ ! -d "$src_dir" ]]; then
        warn "Source directory not found, skipping: $src_dir"
        return 1
    fi
    mkdir -p "$dest_dir"
    cp -a "$src_dir"/. "$dest_dir"/
    return 0
}

# Gzip a file in-place if not already gzipped.
gzip_if_needed() {
    local f="$1"
    if [[ -f "$f" && "${f##*.}" != "gz" ]]; then
        info "Compressing $f"
        gzip "$f"
    fi
}

# Count files staged for a section; record for validation.
record_section() {
    local section="$1" dir="$2"
    local count=0
    if [[ -d "$dir" ]]; then
        count=$(find "$dir" -type f | wc -l)
    fi
    SECTION_FILE_COUNTS+=("${section}:${count}")
    if [[ "$count" -eq 0 ]]; then
        warn "Section '$section' has 0 files staged"
    else
        info "Section '$section': $count files staged"
    fi
}

# ── Begin packaging ─────────────────────────────────────────────────────────
echo "================================================================="
echo " prepare_data_asset.sh  |  Tier $TIER  |  $(date)"
echo " Repo root: $REPO_ROOT"
echo " Output:    $OUTPUT"
echo "================================================================="

mkdir -p "$OUTPUT"

# ── Section 01: Human bulk RNA-seq ──────────────────────────────────────────
echo ""
echo ">>> Section 01: Human bulk RNA-seq"
DEST_01="$OUTPUT/human_bulk_rnaseq"
SRC_01="$REPO_ROOT/data/human_bulk_rnaseq/raw"
safe_copy_contents "$SRC_01" "$DEST_01"
record_section "01_human_bulk_rnaseq" "$DEST_01"

# ── Section 02: Rat snRNA-seq ──────────────────────────────────────────────
echo ""
echo ">>> Section 02: Rat snRNA-seq"
DEST_02="$OUTPUT/rat_snrnaseq"
mkdir -p "$DEST_02"
for i in 1 2 3 4 5 6; do
    nuclei_dir="$REPO_ROOT/data/rat_snrnaseq/raw/Lee_021924_Nuclei${i}"
    if [[ -d "$nuclei_dir" ]]; then
        # Copy the filtered_feature_bc_matrix (barcodes, features, matrix)
        safe_copy "$nuclei_dir" "$DEST_02/Lee_021924_Nuclei${i}" || true
    else
        warn "Nuclei dir not found: $nuclei_dir"
    fi
done
# Include metadata
safe_copy "$REPO_ROOT/data/rat_snrnaseq/metadata/Rat_scRNAseq_AgeGroup.txt" \
          "$DEST_02/Rat_scRNAseq_AgeGroup.txt" || true
record_section "02_rat_snrnaseq" "$DEST_02"

# ── Section 03_comparison: Original vs corrected results ────────────────────
echo ""
echo ">>> Section 03_comparison: Pre/post correction CSVs"
DEST_03="$OUTPUT/03_comparison"
safe_copy_contents "$REPO_ROOT/results/original" "$DEST_03/original" || true
safe_copy_contents "$REPO_ROOT/results/corrected" "$DEST_03/corrected" || true
record_section "03_comparison" "$DEST_03"

# ── Section 03_wes: Rat whole-exome sequencing ─────────────────────────────
echo ""
echo ">>> Section 03_wes: Rat WES data"
DEST_03W="$OUTPUT/rat_wes"
safe_copy_contents "$REPO_ROOT/data/rat_wes/raw" "$DEST_03W" || true
record_section "03_wes" "$DEST_03W"

# ── Section 04: Human scRNA-seq (Xu et al. atlas) ─────────────────────────
echo ""
echo ">>> Section 04: Human scRNA-seq (Xu atlas)"
DEST_04="$OUTPUT/human_scrnaseq"
mkdir -p "$DEST_04"
for f in barcodes.tsv genes.tsv matrix.mtx metadata.csv; do
    safe_copy "$XU_ATLAS/$f" "$DEST_04/$f" || true
done
# Compress matrix.mtx if not already gzipped
if [[ -f "$DEST_04/matrix.mtx" ]]; then
    gzip_if_needed "$DEST_04/matrix.mtx"
fi
record_section "04_human_scrnaseq" "$DEST_04"

# ── Section 05: Rat bulk RNA-seq ───────────────────────────────────────────
echo ""
echo ">>> Section 05: Rat bulk RNA-seq"
DEST_05="$OUTPUT/rat_bulk_rnaseq"
safe_copy_contents "$REPO_ROOT/data/rat_bulk_rnaseq" "$DEST_05" || true
record_section "05_rat_bulk_rnaseq" "$DEST_05"

# ── Section 06: Spatial biopsies (HCC22-088) ───────────────────────────────
echo ""
echo ">>> Section 06: Spatial biopsies"
DEST_06="$OUTPUT/spatial_biopsies"
mkdir -p "$DEST_06"
# The repo has a symlink to the biopsy_adatas.pkl in CITEgeist
if [[ -f "$BIOPSY_PKL" ]]; then
    safe_copy "$BIOPSY_PKL" "$DEST_06/biopsy_adatas.pkl" || true
    # Compress pkl (2.5 GB → ~65 MB gzipped)
    if [[ -f "$DEST_06/biopsy_adatas.pkl" && ! -f "$DEST_06/biopsy_adatas.pkl.gz" ]]; then
        info "Compressing biopsy_adatas.pkl (saves ~2.4 GB)"
        gzip "$DEST_06/biopsy_adatas.pkl"
    fi
else
    # Try the repo symlink
    REPO_SYMLINK="$REPO_ROOT/analysis/06_spatial_biopsies/data/biopsy_adatas.pkl"
    if [[ -L "$REPO_SYMLINK" || -f "$REPO_SYMLINK" ]]; then
        cp -L "$REPO_SYMLINK" "$DEST_06/biopsy_adatas.pkl"
        if [[ -f "$DEST_06/biopsy_adatas.pkl" && ! -f "$DEST_06/biopsy_adatas.pkl.gz" ]]; then
            info "Compressing biopsy_adatas.pkl (saves ~2.4 GB)"
            gzip "$DEST_06/biopsy_adatas.pkl"
        fi
    else
        warn "biopsy_adatas.pkl not found at either location"
    fi
fi
record_section "06_spatial_biopsies" "$DEST_06"

# ── Section 07: Organoid single-cell ───────────────────────────────────────
echo ""
echo ">>> Section 07: Organoid single-cell (Cell Ranger matrices)"
DEST_07="$OUTPUT/organoid_single_cell"
mkdir -p "$DEST_07"
for pool_dir in "$NEIL_PROJECT/pool1_flex/outs/per_sample_outs" \
                "$NEIL_PROJECT/pool2_flex/outs/per_sample_outs"; do
    if [[ ! -d "$pool_dir" ]]; then
        warn "Pool dir not found: $pool_dir"
        continue
    fi
    for sample_dir in "$pool_dir"/OS*; do
        sample_name="$(basename "$sample_dir")"
        matrix_dir="$sample_dir/count/sample_filtered_feature_bc_matrix"
        if [[ -d "$matrix_dir" ]]; then
            safe_copy "$matrix_dir" "$DEST_07/${sample_name}/filtered_feature_bc_matrix" || true
        else
            warn "Filtered matrix not found for $sample_name: $matrix_dir"
        fi
    done
done
record_section "07_organoid_single_cell" "$DEST_07"

# ── GMT gene sets ──────────────────────────────────────────────────────────
echo ""
echo ">>> GMT gene set files"
DEST_GMT="$OUTPUT/gmt"
safe_copy_contents "$REPO_ROOT/data/gmt" "$DEST_GMT" || true
record_section "gmt" "$DEST_GMT"

# ── Tier 2: Precomputed intermediates ──────────────────────────────────────
if [[ "$TIER" == "2" ]]; then
    echo ""
    echo ">>> Tier 2: Precomputed intermediates"
    PRECOMP="$OUTPUT/precomputed"
    mkdir -p "$PRECOMP"

    # 02 outputs: seurat_annotated.rds
    echo "  >> 02: Seurat annotated RDS"
    SEURAT_02="$REPO_ROOT/analysis/02_rat_snrnaseq/outputs/seurat_annotated.rds"
    if [[ -f "$SEURAT_02" ]]; then
        safe_copy "$SEURAT_02" "$PRECOMP/02_outputs/seurat_annotated.rds" || true
    else
        # Copy all outputs from section 02
        safe_copy_contents "$REPO_ROOT/analysis/02_rat_snrnaseq/outputs" "$PRECOMP/02_outputs" || true
    fi
    record_section "precomputed_02" "$PRECOMP/02_outputs"

    # 04 outputs: Seurat objects
    echo "  >> 04: Seurat Xu2024 RDS"
    safe_copy "$REPO_ROOT/data/human_scrnaseq/SeuratObj_Xu2024_HRpos_AfterQCSCT.rds" \
              "$PRECOMP/04_outputs/SeuratObj_Xu2024_HRpos_AfterQCSCT.rds" || true
    record_section "precomputed_04" "$PRECOMP/04_outputs"

    # 07 outputs: h5ad files
    echo "  >> 07: Organoid h5ad files"
    PROCESSED_H5AD="$NEIL_PROJECT/data/processed"
    if [[ -d "$PROCESSED_H5AD" ]]; then
        for h5 in "$PROCESSED_H5AD"/*.h5ad; do
            [[ -f "$h5" ]] || continue
            safe_copy "$h5" "$PRECOMP/07_outputs/$(basename "$h5")" || true
        done
    else
        warn "Processed h5ad directory not found: $PROCESSED_H5AD"
    fi
    record_section "precomputed_07" "$PRECOMP/07_outputs"
fi

# ── Generate data_README.md with SHA256 checksums ──────────────────────────
echo ""
echo ">>> Generating data_README.md with checksums..."
README="$OUTPUT/data_README.md"
{
    echo "# ERpBRCA_OlderWomen Data Asset"
    echo ""
    echo "Generated: $(date -Iseconds)"
    echo "Tier: $TIER"
    echo "Source repo: $REPO_ROOT"
    echo ""
    echo "## Section file counts"
    echo ""
    for entry in "${SECTION_FILE_COUNTS[@]}"; do
        section="${entry%%:*}"
        count="${entry##*:}"
        echo "- **${section}**: ${count} files"
    done
    echo ""
    echo "## SHA256 Checksums"
    echo ""
    echo '```'
} > "$README"

# Compute checksums (sorted for reproducibility)
find "$OUTPUT" -type f ! -name "data_README.md" -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    | sed "s|${OUTPUT}/||" \
    >> "$README"

echo '```' >> "$README"

# ── Validation summary ─────────────────────────────────────────────────────
echo ""
echo "================================================================="
echo " Validation Summary"
echo "================================================================="
EMPTY_SECTIONS=0
for entry in "${SECTION_FILE_COUNTS[@]}"; do
    section="${entry%%:*}"
    count="${entry##*:}"
    if [[ "$count" -eq 0 ]]; then
        echo "  [FAIL] $section: 0 files"
        ((EMPTY_SECTIONS++)) || true
    else
        echo "  [ OK ] $section: $count files"
    fi
done

# ── Total size ──────────────────────────────────────────────────────────────
echo ""
TOTAL_SIZE=$(du -sh "$OUTPUT" | cut -f1)
echo "Total staging size: $TOTAL_SIZE"
echo "Warnings: $WARNINGS"
echo "Empty sections: $EMPTY_SECTIONS"
echo ""

if [[ "$EMPTY_SECTIONS" -gt 0 ]]; then
    echo "[WARNING] $EMPTY_SECTIONS section(s) have zero files. Review warnings above."
fi

echo "Staging directory: $OUTPUT"
echo "Done."
