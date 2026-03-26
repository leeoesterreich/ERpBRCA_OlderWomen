#!/bin/bash
#SBATCH --job-name=cr_pool2
#SBATCH --output=logs/cellranger_pool2_%j.out
#SBATCH --error=logs/cellranger_pool2_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Cell Ranger multi for Pool 2 (OS05-07: E2, E2+ICI, E2+HSD17B7i)
# 10X Flex demultiplexing - each pool must run separately

module load cellranger/9.0.1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Starting Cell Ranger multi POOL 2 at $(date)"
echo "Working directory: $(pwd)"
echo "Samples: OS05 (E2), OS06 (E2+ICI), OS07 (E2+HSD17B7i)"

cellranger multi \
    --id=pool2_flex \
    --csv=configs/cellranger_multi_config_pool2.csv \
    --localcores=16 \
    --localmem=120

echo "Cell Ranger multi POOL 2 completed at $(date)"
