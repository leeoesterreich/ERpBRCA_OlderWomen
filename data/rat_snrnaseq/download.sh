#!/bin/bash
# data/rat_snrnaseq/download.sh
# Downloads GSE276758 rat snRNA-seq data
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"
META_DIR="$SCRIPT_DIR/metadata"

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
wget -q -P "$RAW_DIR" "${GEO_BASE}/GSE276758_RAW.tar" || echo "Warning: Could not download 10X data"

# Extract if downloaded
if [[ -f "$RAW_DIR/GSE276758_RAW.tar" ]]; then
    tar -xf "$RAW_DIR/GSE276758_RAW.tar" -C "$RAW_DIR"
    rm "$RAW_DIR/GSE276758_RAW.tar"
fi

# Create age group annotation file
cat > "$META_DIR/Rat_scRNAseq_AgeGroup.txt" << 'EOF'
CaseID	AgeGroup
Lee_021924_Nuclei1	Aged
Lee_021924_Nuclei2	Aged
Lee_021924_Nuclei3	Aged
Lee_021924_Nuclei4	Young
Lee_021924_Nuclei5	Young
Lee_021924_Nuclei6	Young
EOF

echo "Rat snRNA-seq download complete."
