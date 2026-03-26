#!/bin/bash
# data/rat_snrnaseq/download.sh
# Downloads GSE276758 rat snRNA-seq data
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"
META_DIR="$SCRIPT_DIR/metadata"

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

# Helper function to verify directory exists and is non-empty
verify_dir() {
    local dir="$1"
    local desc="${2:-directory}"
    if [[ ! -d "$dir" ]]; then
        echo "ERROR: $desc not found: $dir"
        exit 1
    fi
    if [[ -z "$(ls -A "$dir")" ]]; then
        echo "ERROR: $desc is empty: $dir"
        exit 1
    fi
    echo "  Verified directory: $dir"
}

mkdir -p "$RAW_DIR" "$META_DIR"

echo "Downloading GSE276758 data..."

# Check if files already exist
if [[ -d "$RAW_DIR/Lee_021924_Nuclei1" ]]; then
    echo "Data already downloaded, skipping..."
    exit 0
fi

# Note: 10X data from GEO requires specific download structure
# Update with actual GEO supplementary file URLs
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE276nnn/GSE276758/suppl"

# Download the tar archive of 10X matrices
echo "  Downloading 10X data archive..."
wget --show-progress -P "$RAW_DIR" "${GEO_BASE}/GSE276758_RAW.tar"
verify_file "$RAW_DIR/GSE276758_RAW.tar" "10X data archive"

# Extract the archive
echo "  Extracting archive..."
tar -xvf "$RAW_DIR/GSE276758_RAW.tar" -C "$RAW_DIR"

# Verify extraction produced expected content
# Check for at least one expected directory
if [[ ! -d "$RAW_DIR/Lee_021924_Nuclei1" ]]; then
    # Some GEO archives extract files directly, check for any .gz files
    if ! ls "$RAW_DIR"/*.gz 1>/dev/null 2>&1; then
        echo "ERROR: Extraction did not produce expected files"
        exit 1
    fi
    echo "  Verified: Archive extracted files successfully"
else
    verify_dir "$RAW_DIR/Lee_021924_Nuclei1" "extracted 10X data"
fi

# Only remove archive after successful extraction verification
echo "  Removing archive..."
rm "$RAW_DIR/GSE276758_RAW.tar"

# Create age group annotation file
echo "  Creating metadata file..."
cat > "$META_DIR/Rat_scRNAseq_AgeGroup.txt" << 'EOF'
CaseID	AgeGroup
Lee_021924_Nuclei1	Aged
Lee_021924_Nuclei2	Aged
Lee_021924_Nuclei3	Aged
Lee_021924_Nuclei4	Young
Lee_021924_Nuclei5	Young
Lee_021924_Nuclei6	Young
EOF

verify_file "$META_DIR/Rat_scRNAseq_AgeGroup.txt" "metadata file"

echo "Rat snRNA-seq download complete."
