#!/bin/bash
# Convert all methods markdown files to Word docx format
# Requires: pandoc (available via module load on CRC)
# Usage: bash scripts/methods_to_docx.sh

set -eo pipefail

METHODS_DIR="docs/methods"
OUTPUT_DIR="docs/methods/docx"
mkdir -p "$OUTPUT_DIR"

# Load pandoc if available as module
module load pandoc 2>/dev/null || true

# Check pandoc is available
if ! command -v pandoc &>/dev/null; then
    echo "ERROR: pandoc not found. Try: module load pandoc"
    exit 1
fi

echo "Converting methods .md files to .docx..."

for md_file in "$METHODS_DIR"/*.md; do
    basename=$(basename "$md_file" .md)
    # Skip non-methods files
    [[ "$basename" == "methods_vs_code_audit" ]] && continue

    docx_file="$OUTPUT_DIR/${basename}.docx"
    echo "  $basename.md -> $basename.docx"
    pandoc "$md_file" \
        -f markdown \
        -t docx \
        --reference-doc="$OUTPUT_DIR/reference.docx" 2>/dev/null \
        -o "$docx_file" \
        || pandoc "$md_file" -f markdown -t docx -o "$docx_file"
done

echo "Done. Output in $OUTPUT_DIR/"
ls -la "$OUTPUT_DIR"/*.docx 2>/dev/null

# Also create a single combined document
echo ""
echo "Creating combined methods document..."
COMBINED="$OUTPUT_DIR/all_methods_combined.docx"
# Concatenate all md files in order, with section breaks
cat "$METHODS_DIR"/00_main_text_methods.md \
    "$METHODS_DIR"/00b_supplementary_methods.md \
    "$METHODS_DIR"/01_human_bulk_rnaseq_methods.md \
    "$METHODS_DIR"/02_rat_snrnaseq_methods.md \
    "$METHODS_DIR"/03_rat_wes_methods.md \
    "$METHODS_DIR"/04_human_scrnaseq_methods.md \
    "$METHODS_DIR"/05_rat_bulk_rnaseq_methods.md \
    "$METHODS_DIR"/06_spatial_biopsies_methods.md \
    "$METHODS_DIR"/07_organoid_scrnaseq_methods.md \
    2>/dev/null | pandoc -f markdown -t docx -o "$COMBINED"
echo "  Combined: $COMBINED"
