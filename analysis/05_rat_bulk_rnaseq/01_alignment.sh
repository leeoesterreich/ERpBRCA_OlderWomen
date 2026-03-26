#!/bin/bash
#SBATCH --job-name=rat_alignment
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00
#SBATCH --mem=128G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/alignment_%j.out
#SBATCH --error=logs/alignment_%j.err

# Rat bulk RNA-seq alignment with STAR
# Inputs: Trimmed FASTQ files from fastp
# Outputs: Sorted BAM files

set -euo pipefail

module purge
module load gcc/8.2.0
module load star/2.7.11b

# Paths - use absolute paths for SLURM compatibility
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/outputs/trimmed"
OUTPUT_DIR="${SCRIPT_DIR}/outputs/aligned"
# Reference genome directory (set to your local path)
# Index built with STAR v2.7.4a; binary is v2.7.11b (minor version mismatch, backward compatible)
GENOME_DIR="${GENOME_DIR:-data/external/rat_genome}"
THREADS=64

mkdir -p "$OUTPUT_DIR"
mkdir -p "${SCRIPT_DIR}/logs"

echo "=== STAR Alignment ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Genome: $GENOME_DIR"

# Verify genome index exists
if [ ! -f "${GENOME_DIR}/SAindex" ]; then
    echo "ERROR: STAR index not found at $GENOME_DIR"
    exit 1
fi

# Align each sample
for read1 in "$INPUT_DIR"/*_R1_trimmed.fastq.gz; do
    [ -f "$read1" ] || continue

    sample=$(basename "$read1" _R1_trimmed.fastq.gz)
    read2="${INPUT_DIR}/${sample}_R2_trimmed.fastq.gz"
    output_prefix="${OUTPUT_DIR}/${sample}_"

    if [ ! -f "$read2" ]; then
        echo "Warning: R2 not found for $sample, skipping"
        continue
    fi

    # Skip if already aligned
    if [ -f "${output_prefix}Aligned.sortedByCoord.out.bam" ]; then
        echo "Skipping $sample (already aligned)"
        continue
    fi

    echo "Aligning: $sample"
    STAR --genomeDir "$GENOME_DIR" \
         --readFilesIn "$read1" "$read2" \
         --readFilesCommand zcat \
         --runThreadN $THREADS \
         --outFileNamePrefix "$output_prefix" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts \
         --outBAMsortingThreadN 8

    echo "Completed: $sample"
done

echo "=== Alignment complete ==="
echo "Output files:"
ls -la "$OUTPUT_DIR"/*.bam 2>/dev/null || echo "No BAM files yet"
