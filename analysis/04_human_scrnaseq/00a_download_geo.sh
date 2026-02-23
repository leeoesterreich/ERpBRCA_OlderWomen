#!/bin/bash
#SBATCH --job-name=download_geo
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err

# 00a_download_geo.sh
# Download Wu et al. (GSE176078) scRNA-seq count matrices from GEO
#
# Outputs:
#   - data/human_scrnaseq/raw/CID*/matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data/human_scrnaseq/raw"

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

echo "=== Downloading Wu et al. GSE176078 scRNA-seq data ==="
echo "Target directory: $DATA_DIR"

# Download main tar file from GEO
GEO_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
TAR_FILE="GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"

if [[ ! -f "$TAR_FILE" ]]; then
    echo "Downloading $TAR_FILE..."
    wget -c "$GEO_URL"
else
    echo "$TAR_FILE already exists, skipping download"
fi

# Extract
echo "Extracting archive..."
tar -xzf "$TAR_FILE"

# List extracted contents
echo "Extracted files:"
ls -la

# Verify expected sample directories exist
EXPECTED_SAMPLES="CID3941 CID3948 CID4040 CID4067 CID4290A CID4461 CID4463 CID4471 CID4530N CID4535"
for sample in $EXPECTED_SAMPLES; do
    if [[ -d "$sample" ]] || [[ -d "$DATA_DIR/$sample" ]]; then
        echo "  Found: $sample"
    else
        echo "  WARNING: Missing sample directory: $sample"
    fi
done

echo "=== Download complete ==="
