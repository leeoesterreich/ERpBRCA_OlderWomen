#!/bin/bash
#SBATCH --job-name=build_cellranger_ref
#SBATCH --output=logs/build_ref_%j.out
#SBATCH --error=logs/build_ref_%j.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Build Cell Ranger reference GRCh38-2024-A
# Based on 10x Genomics official build script

module load cellranger/9.0.1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Genome metadata
genome="GRCh38"
version="2024-A"

# Set up source and build directories
build="references/GRCh38-GENCODEv44_build"
mkdir -p "$build"

# Download source files if they do not exist in reference_sources/ folder
source="references/reference_sources"
mkdir -p "$source"

# Using release 109 FASTA for GRCh38 instead of release 110 FASTA -- release 110 moved from GRCh38.p13 to GRCh38.p14,
# which unmasked the pseudo-autosomal region. This causes ambiguous mappings to PAR locus genes.
# No other sequence changes were made to the primary assembly.
fasta_url="http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.v44.primary_assembly.annotation.gtf"

echo "Downloading FASTA..."
if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi

echo "Downloading GTF..."
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi

# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
echo "Modifying FASTA headers..."
fasta_modified="$build/$(basename "$fasta_in").modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"

# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
echo "Modifying GTF..."
gtf_modified="$build/$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

# Define string patterns for GTF tags
# NOTES:
# Since Ensembl 110, polymorphic pseudogenes are now just protein_coding.
# Readthrough genes are annotated with the readthrough_transcript tag.
BIOTYPE_PATTERN=\
"(protein_coding|protein_coding_LoF|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""

# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
echo "Building gene allowlist..."
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"

# NOTES:
# Since Ensembl 110, the PAR locus genes are included on chrY as copies of chrX
# Using the GRCh38.p13 assembly hard masks these regions on chrY, but removing the
# chrY PAR genes is still desirable so they do not end up as extra entries in the output.
# The awk command below excludes all PAR_Y genes, including XGY2.
# The non-coding gene XGY2 straddles the PAR1 boundary on chrY, and is homologous to XG on chrX.
# GRCh38-2024-A excludes XGY2, but includes SRY and ENSG00000286130, which are in an intron of XGY2,
# and RPS4Y1, which overlaps XGY2.

# Filter the GTF file based on the gene allowlist
echo "Filtering GTF..."
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
# Filter to the gene allowlist, and then remove PAR_Y genes
grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    | awk -F "\t" '$1 != "chrY" || $1 == "chrY" && $4 >= 2752083 && $4 < 56887903 && !/ENSG00000290840/' \
    >> "$gtf_filtered"

# Create reference package
echo "Running cellranger mkref..."
cellranger mkref --ref-version="$version" \
    --genome="$genome" --fasta="$fasta_modified" --genes="$gtf_filtered" \
    --nthreads=16

# Move output to references folder
mv GRCh38 references/refdata-gex-GRCh38-2024-A-custom

echo "Reference build complete!"
echo "Output: references/refdata-gex-GRCh38-2024-A-custom"
