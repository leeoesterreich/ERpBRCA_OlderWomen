#!/bin/bash
# Run human scRNA-seq analysis pipeline
# Usage: bash run_analysis.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Human scRNA-seq Pipeline ==="
echo "Directory: $SCRIPT_DIR"

# Submit SLURM job
sbatch run_analysis.sbatch

echo "=== Job submitted ==="
