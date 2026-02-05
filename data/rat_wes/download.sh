#!/bin/bash
# data/rat_wes/download.sh
# Downloads GSE276759 rat WES data
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"

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

echo "  Downloading WES data archive..."
wget --show-progress -P "$RAW_DIR" "${GEO_BASE}/GSE276759_processed_data.tar.gz"
verify_file "$RAW_DIR/GSE276759_processed_data.tar.gz" "WES data archive"

# Count files before extraction to verify it worked
echo "  Extracting archive..."
files_before=$(find "$RAW_DIR" -type f | wc -l)
tar -xzvf "$RAW_DIR/GSE276759_processed_data.tar.gz" -C "$RAW_DIR"
files_after=$(find "$RAW_DIR" -type f | wc -l)

# Verify extraction produced new files
if [[ "$files_after" -le "$files_before" ]]; then
    echo "ERROR: Extraction did not produce any new files"
    exit 1
fi
echo "  Verified: Extracted $((files_after - files_before)) new files"

# Only remove archive after successful extraction verification
echo "  Removing archive..."
rm "$RAW_DIR/GSE276759_processed_data.tar.gz"

# Mark as downloaded only after all steps succeed
touch "$RAW_DIR/.downloaded"
echo "Rat WES download complete."
