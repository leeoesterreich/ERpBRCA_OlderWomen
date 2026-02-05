#!/bin/bash
# data/download_all.sh
# Master script to download all data from GEO
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=== Downloading all GEO data ==="

echo "[1/3] Downloading human bulk RNA-seq (GSE276755)..."
bash "$SCRIPT_DIR/human_bulk_rnaseq/download.sh"

echo "[2/3] Downloading rat snRNA-seq (GSE276758)..."
bash "$SCRIPT_DIR/rat_snrnaseq/download.sh"

echo "[3/3] Downloading rat WES (GSE276759)..."
bash "$SCRIPT_DIR/rat_wes/download.sh"

echo "=== All downloads complete ==="
