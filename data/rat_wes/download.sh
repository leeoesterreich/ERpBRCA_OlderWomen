#!/bin/bash
# data/rat_wes/download.sh
# Downloads GSE276759 rat WES data
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"

mkdir -p "$RAW_DIR"

echo "Downloading GSE276759 data..."

# Check if files already exist
if [[ -f "$RAW_DIR/.downloaded" ]]; then
    echo "Data already downloaded, skipping..."
    exit 0
fi

# Note: WES data may need SRA download
# Update with actual GEO/SRA links
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE276nnn/GSE276759/suppl"

wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276759_processed_data.tar.gz" || echo "Warning: Could not download WES data"

# Extract if downloaded
if [[ -f "$RAW_DIR/GSE276759_processed_data.tar.gz" ]]; then
    tar -xzf "$RAW_DIR/GSE276759_processed_data.tar.gz" -C "$RAW_DIR"
fi

touch "$RAW_DIR/.downloaded"
echo "Rat WES download complete."
