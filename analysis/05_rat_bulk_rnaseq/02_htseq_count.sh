#!/bin/bash
#SBATCH --job-name=rat_htseq
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/htseq_%j.out
#SBATCH --error=logs/htseq_%j.err

# HTSeq read counting for rat bulk RNA-seq
# Inputs: Aligned BAM files
# Outputs: Gene count files

set -euo pipefail

module purge
module load htseq/0.13.5

# Paths - use absolute paths for SLURM compatibility
SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/05_rat_bulk_rnaseq"
INPUT_DIR="${SCRIPT_DIR}/outputs/aligned"
OUTPUT_DIR="${SCRIPT_DIR}/outputs/counts"
GTF_FILE="/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Rat/Rattus_norvegicus.mRatBN7.2.112.gtf"

# Strand setting - reverse for typical Illumina TruSeq
STRAND="reverse"

mkdir -p "$OUTPUT_DIR"

echo "=== HTSeq Count ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "GTF: $GTF_FILE"
echo "Strand: $STRAND"

for bam_file in "$INPUT_DIR"/*Aligned.sortedByCoord.out.bam; do
    [ -f "$bam_file" ] || continue

    sample=$(basename "$bam_file" _Aligned.sortedByCoord.out.bam)
    output_file="${OUTPUT_DIR}/${sample}.txt"

    # Skip if already counted
    if [ -f "$output_file" ]; then
        echo "Skipping $sample (already counted)"
        continue
    fi

    echo "Counting: $sample"
    htseq-count \
        -f bam \
        -r pos \
        -s "$STRAND" \
        -t exon \
        -i gene_id \
        "$bam_file" \
        "$GTF_FILE" > "$output_file"

    echo "Created: $output_file"
done

echo "=== HTSeq complete ==="
echo "Output files:"
ls -la "$OUTPUT_DIR"
