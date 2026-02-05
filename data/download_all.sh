#!/bin/bash
# data/download_all.sh
# Master script to download all data from GEO
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Downloading all GEO data ==="

echo "[1/3] Downloading human bulk RNA-seq (GSE276755)..."
if ! bash "$SCRIPT_DIR/human_bulk_rnaseq/download.sh"; then
    echo "ERROR: Human bulk RNA-seq download failed"
    exit 1
fi

echo "[2/3] Downloading rat snRNA-seq (GSE276758)..."
if ! bash "$SCRIPT_DIR/rat_snrnaseq/download.sh"; then
    echo "ERROR: Rat snRNA-seq download failed"
    exit 1
fi

echo "[3/3] Downloading rat WES (GSE276759)..."
if ! bash "$SCRIPT_DIR/rat_wes/download.sh"; then
    echo "ERROR: Rat WES download failed"
    exit 1
fi

echo "=== All downloads complete ==="
