#!/bin/bash
# data/human_bulk_rnaseq/download.sh
# Downloads GSE276755 human bulk RNA-seq data
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"
EXT_DIR="$SCRIPT_DIR/external"

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

wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276755_HumanERpAge_39404g168s_FeatureCount.txt.gz" || echo "Warning: Could not download count file"
wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276755_HumanERpAge_39404g168s_TPMlog2.txt.gz" || echo "Warning: Could not download TPM file"
wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276755_HumanERpAge_BulkRNAseq_SampleInformation.txt.gz" || echo "Warning: Could not download sample info"

# Decompress if downloaded
gunzip -k "$RAW_DIR"/*.gz 2>/dev/null || true

# Download gene annotation from NCBI
echo "Downloading gene annotation..."
wget -q -P "$EXT_DIR" "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
gunzip -k "$EXT_DIR/Homo_sapiens.gene_info.gz" 2>/dev/null || true
mv "$EXT_DIR/Homo_sapiens.gene_info" "$EXT_DIR/Homo_sapiens.gene_info.txt" 2>/dev/null || true

echo "Human bulk RNA-seq download complete."
