#!/bin/bash
#SBATCH --job-name=rat_htseq
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# HTSeq read counting for rat bulk RNA-seq
# Inputs: Aligned BAM files
# Outputs: Gene count files

set -euo pipefail

module purge
module load htseq/0.13.5

PROJECT_ROOT="${PROJECT_ROOT:-$(dirname $(dirname $(realpath $0)))}"
INPUT_DIR="${INPUT_DIR:-$PROJECT_ROOT/analysis/05_rat_bulk_rnaseq/outputs/aligned}"
OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT_ROOT/analysis/05_rat_bulk_rnaseq/outputs/counts}"
GTF_FILE="${GTF_FILE:-/path/to/Rattus_norvegicus.mRatBN7.2.112.gtf}"

# FIX: Added strand specification - check library prep protocol
# Options: yes (forward), reverse, no (unstranded)
STRAND="${STRAND:-reverse}"

mkdir -p "$OUTPUT_DIR"

echo "=== HTSeq Count ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Strand: $STRAND"

for bam_file in "$INPUT_DIR"/*Aligned.sortedByCoord.out.bam; do
    [ -f "$bam_file" ] || continue

    sample=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)
    output_file="${OUTPUT_DIR}/${sample}.txt"

    echo "Counting: $sample"
    # FIX: Added -s flag for strand specification
    htseq-count -f bam -r pos -s "$STRAND" "$bam_file" "$GTF_FILE" > "$output_file"

    echo "Created: $output_file"
done

echo "=== HTSeq complete ==="
