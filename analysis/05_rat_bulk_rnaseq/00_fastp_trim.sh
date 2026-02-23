#!/bin/bash
#SBATCH --job-name=rat_fastp
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 4:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/fastp_%j.out
#SBATCH --error=logs/fastp_%j.err

# Rat bulk RNA-seq trimming with fastp
# Inputs: Raw FASTQ files
# Outputs: Trimmed FASTQ files + QC reports

set -euo pipefail

module purge
module load fastp/0.23.4

# Paths - use absolute paths for SLURM compatibility
RAW_DIR="/ix1/alee/LO_LAB/General/Lab_Data/20240622_Neil_Rahul_RNASeqRat"
SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/05_rat_bulk_rnaseq"
OUTPUT_DIR="${SCRIPT_DIR}/outputs/trimmed"

mkdir -p "$OUTPUT_DIR"
mkdir -p "${SCRIPT_DIR}/logs"

echo "=== fastp Trimming ==="
echo "Input: $RAW_DIR"
echo "Output: $OUTPUT_DIR"

# Process each sample
for read1 in "$RAW_DIR"/*_R1_001.fastq.gz; do
    [ -f "$read1" ] || continue

    # Extract sample name (e.g., 102-FF-Tumor from 102-FF-Tumor_R1_001.fastq.gz)
    sample=$(basename "$read1" _R1_001.fastq.gz)
    read2="${RAW_DIR}/${sample}_R2_001.fastq.gz"

    if [ ! -f "$read2" ]; then
        echo "Warning: R2 not found for $sample, skipping"
        continue
    fi

    echo "Processing: $sample"

    fastp \
        -i "$read1" \
        -I "$read2" \
        -o "${OUTPUT_DIR}/${sample}_R1_trimmed.fastq.gz" \
        -O "${OUTPUT_DIR}/${sample}_R2_trimmed.fastq.gz" \
        --html "${OUTPUT_DIR}/${sample}_fastp.html" \
        --json "${OUTPUT_DIR}/${sample}_fastp.json" \
        --thread 8 \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 36

    echo "Completed: $sample"
done

echo "=== Trimming complete ==="
echo "Output files:"
ls -la "$OUTPUT_DIR"
