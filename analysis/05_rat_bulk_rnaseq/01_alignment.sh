#!/bin/bash
#SBATCH --job-name=rat_alignment
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -t 1-00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Rat bulk RNA-seq alignment with STAR
# Inputs: Trimmed FASTQ files
# Outputs: Sorted BAM files

set -euo pipefail

module purge
module load gcc/8.2.0
module load star/2.7.9a

# Configurable paths
PROJECT_ROOT="${PROJECT_ROOT:-$(dirname $(dirname $(realpath $0)))}"
INPUT_DIR="${INPUT_DIR:-$PROJECT_ROOT/data/rat_bulk_rnaseq/trimmed}"
OUTPUT_DIR="${OUTPUT_DIR:-$PROJECT_ROOT/analysis/05_rat_bulk_rnaseq/outputs/aligned}"
GENOME_DIR="${GENOME_DIR:-/path/to/Reference_genome/Rat}"
GENOME_FASTA="${GENOME_DIR}/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa"
GTF_FILE="${GENOME_DIR}/Rattus_norvegicus.mRatBN7.2.112.gtf"
THREADS=64

mkdir -p "$OUTPUT_DIR"

echo "=== STAR Alignment ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"

# Genome indexing (skip if index exists)
if [ ! -f "${GENOME_DIR}/SAindex" ]; then
    echo "Indexing reference genome..."
    STAR --runMode genomeGenerate \
         --genomeDir "$GENOME_DIR" \
         --genomeFastaFiles "$GENOME_FASTA" \
         --sjdbGTFfile "$GTF_FILE" \
         --runThreadN $THREADS
fi

# Align each sample
# FIX: Changed input_dir1 to INPUT_DIR
for read1 in "$INPUT_DIR"/*_R1_paired.fastq; do
    [ -f "$read1" ] || continue

    sample=$(basename "$read1" _R1_paired.fastq)
    read2="$INPUT_DIR/${sample}_R2_paired.fastq"
    output_prefix="$OUTPUT_DIR/${sample}_"

    if [ ! -f "$read2" ]; then
        echo "Warning: R2 not found for $sample, skipping"
        continue
    fi

    echo "Aligning: $sample"
    STAR --genomeDir "$GENOME_DIR" \
         --readFilesIn "$read1" "$read2" \
         --runThreadN $THREADS \
         --outFileNamePrefix "$output_prefix" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts

    echo "Completed: $sample"
done

echo "=== Alignment complete ==="
