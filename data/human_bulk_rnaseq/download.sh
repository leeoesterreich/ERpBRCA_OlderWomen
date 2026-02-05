#!/bin/bash
# data/human_bulk_rnaseq/download.sh
# Downloads GSE276755 human bulk RNA-seq data
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"
EXT_DIR="$SCRIPT_DIR/external"

# Helper function to verify file exists and is non-empty
verify_file() {
    local file="$1"
    local desc="${2:-file}"
    if [[ ! -f "$file" ]]; then
        echo "ERROR: $desc not found: $file"
        exit 1
    fi
    if [[ ! -s "$file" ]]; then
        echo "ERROR: $desc is empty: $file"
        exit 1
    fi
    echo "  Verified: $file"
}

mkdir -p "$RAW_DIR" "$EXT_DIR"

echo "Downloading GSE276755 data..."

# Check if files already exist
if [[ -f "$RAW_DIR/HumanERpAge_39404g168s_FeatureCount.txt" ]]; then
    echo "Data already downloaded, skipping..."
    exit 0
fi

# Download from GEO supplementary files
# Note: Update these URLs with actual GEO links when available
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE276nnn/GSE276755/suppl"

echo "  Downloading count file..."
wget --show-progress -P "$RAW_DIR" "${GEO_BASE}/GSE276755_HumanERpAge_39404g168s_FeatureCount.txt.gz"
verify_file "$RAW_DIR/GSE276755_HumanERpAge_39404g168s_FeatureCount.txt.gz" "count file"

echo "  Downloading TPM file..."
wget --show-progress -P "$RAW_DIR" "${GEO_BASE}/GSE276755_HumanERpAge_39404g168s_TPMlog2.txt.gz"
verify_file "$RAW_DIR/GSE276755_HumanERpAge_39404g168s_TPMlog2.txt.gz" "TPM file"

echo "  Downloading sample info..."
wget --show-progress -P "$RAW_DIR" "${GEO_BASE}/GSE276755_HumanERpAge_BulkRNAseq_SampleInformation.txt.gz"
verify_file "$RAW_DIR/GSE276755_HumanERpAge_BulkRNAseq_SampleInformation.txt.gz" "sample info"

# Decompress downloaded files
echo "  Decompressing files..."
for gz_file in "$RAW_DIR"/*.gz; do
    if [[ -f "$gz_file" ]]; then
        gunzip -k "$gz_file"
        # Verify decompressed file
        decompressed="${gz_file%.gz}"
        verify_file "$decompressed" "decompressed $(basename "$decompressed")"
    fi
done

# Download gene annotation from NCBI
echo "Downloading gene annotation..."
wget --show-progress -P "$EXT_DIR" "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
verify_file "$EXT_DIR/Homo_sapiens.gene_info.gz" "gene annotation archive"

gunzip -k "$EXT_DIR/Homo_sapiens.gene_info.gz"
verify_file "$EXT_DIR/Homo_sapiens.gene_info" "decompressed gene annotation"

mv "$EXT_DIR/Homo_sapiens.gene_info" "$EXT_DIR/Homo_sapiens.gene_info.txt"
verify_file "$EXT_DIR/Homo_sapiens.gene_info.txt" "renamed gene annotation"

echo "Human bulk RNA-seq download complete."
