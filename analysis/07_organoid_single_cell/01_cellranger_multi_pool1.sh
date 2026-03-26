#!/bin/bash
#SBATCH --job-name=cr_pool1
#SBATCH --output=logs/cellranger_pool1_%j.out
#SBATCH --error=logs/cellranger_pool1_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Cell Ranger multi for Pool 1 (OS01-04: Vehicle, E1, E1+ICI, E1+HSD17B7i)
# 10X Flex demultiplexing - each pool must run separately

module load cellranger/9.0.1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Starting Cell Ranger multi POOL 1 at $(date)"
echo "Working directory: $(pwd)"
echo "Samples: OS01 (Vehicle), OS02 (E1), OS03 (E1+ICI), OS04 (E1+HSD17B7i)"

cellranger multi \
    --id=pool1_flex \
    --csv=configs/cellranger_multi_config_pool1.csv \
    --localcores=16 \
    --localmem=120

echo "Cell Ranger multi POOL 1 completed at $(date)"
