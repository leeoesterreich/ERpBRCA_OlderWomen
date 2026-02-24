#!/bin/bash
#SBATCH --job-name=vep_rat
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=logs/vep_%j.out
#SBATCH --error=logs/vep_%j.err

# VEP Annotation for Rat WES VCF files
# Annotates Mutect2 VCF files with variant consequences

set -eo pipefail

module load ensembl-vep/114.2

SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/03_rat_wes"
VCF_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/RatWES_Mutect2_VCF_Input"
VEP_CACHE="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/rattus_norvegicus"
OUTPUT_DIR="$SCRIPT_DIR/outputs/vep"

cd "$SCRIPT_DIR"
mkdir -p logs "$OUTPUT_DIR"

echo "=== Running VEP Annotation ==="
echo "VCF directory: $VCF_DIR"
echo "VEP cache: $VEP_CACHE"
echo "Output: $OUTPUT_DIR"

# Process each VCF file
for vcf in "$VCF_DIR"/*.vcf; do
    sample=$(basename "$vcf" .vcf)

    # Convert sample name to match expected format
    # Rat_O_102 -> 102_tumor_vs_102_spleen
    sample_num=$(echo "$sample" | grep -oE '[0-9]+')
    output_name="${sample_num}_tumor_vs_${sample_num}_spleen"

    echo "Processing: $sample -> $output_name"

    # Create output directory structure matching original
    mkdir -p "$OUTPUT_DIR/$output_name"

    vep \
        --input_file "$vcf" \
        --output_file "$OUTPUT_DIR/$output_name/${output_name}.mutect2.txt" \
        --cache \
        --dir_cache "$VEP_CACHE" \
        --species rattus_norvegicus \
        --assembly Rnor_6.0 \
        --cache_version 95 \
        --offline \
        --force_overwrite \
        --tab \
        --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,Extra" \
        --fork 4

    echo "Completed: $output_name"
done

echo "=== VEP Annotation Complete ==="
echo "Outputs in: $OUTPUT_DIR"
