#!/usr/bin/env bash
# adapt_capsule_paths.sh — Apply CodeOcean path transformations to capsule scripts
#
# Run AFTER sync_to_capsule.sh to replace HPC paths with capsule-relative paths.
# This script is idempotent — safe to run multiple times.
#
# Usage:
#   ./scripts/adapt_capsule_paths.sh                    # Apply to default capsule dir
#   ./scripts/adapt_capsule_paths.sh --target /path      # Override capsule dir
#   ./scripts/adapt_capsule_paths.sh --dry-run           # Show what would change

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CAPSULE_DIR="$(cd "$REPO_ROOT/../ERpBRCA_CodeOcean" 2>/dev/null && pwd)"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --target) CAPSULE_DIR="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

CODE="$CAPSULE_DIR/code"
if [ ! -d "$CODE" ]; then
  echo "ERROR: Capsule code dir not found: $CODE"
  exit 1
fi

CHANGED=0

do_sed() {
  local file="$1" pattern="$2" replacement="$3"
  if grep -q "$pattern" "$file" 2>/dev/null; then
    if [ "$DRY_RUN" = "1" ]; then
      echo "  WOULD FIX: $file"
      grep -n "$pattern" "$file" | head -3
    else
      sed -i "s|$pattern|$replacement|g" "$file"
      CHANGED=$((CHANGED + 1))
    fi
  fi
}

echo "=== CodeOcean Path Adaptation ==="
echo "Capsule: $CAPSULE_DIR"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# 1. SOURCE paths.R at the top of R scripts that use project_root or data paths
# ─────────────────────────────────────────────────────────────────────────────
echo "--- Step 1: Inject paths.R sourcing into R scripts ---"

for rscript in $(find "$CODE" -name "*.R" -not -path "*/config/*" -not -path "*/lib/*"); do
  # Skip if already sources paths.R
  if grep -q 'source.*config.*paths\.R' "$rscript" 2>/dev/null; then
    continue
  fi
  # Only inject if script uses project_root, data_dir, output_dir, or hardcoded /ix1/
  if grep -qE 'project_root|data_dir.*file\.path|output_dir.*file\.path\(script_dir|/ix1/' "$rscript" 2>/dev/null; then
    if [ "$DRY_RUN" = "1" ]; then
      echo "  WOULD INJECT paths.R source: $(basename "$rscript")"
    else
      # Determine depth: mica/ scripts need ../../config, others need ../config
      if echo "$rscript" | grep -q "/mica/"; then
        CONFIG_REL='file.path(script_dir, "..", "..", "config", "paths.R")'
      else
        CONFIG_REL='file.path(script_dir, "..", "config", "paths.R")'
      fi
      # Insert after the script_dir assignment
      sed -i "/^script_dir <- get_script_dir()/a\\
source($CONFIG_REL)" "$rscript"
      CHANGED=$((CHANGED + 1))
      echo "  Injected paths.R: $(basename "$rscript")"
    fi
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
# 2. Replace project_root with ROOT_DIR in R scripts
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- Step 2: Replace project_root → ROOT_DIR ---"

for rscript in $(grep -rl 'project_root' "$CODE" --include="*.R" 2>/dev/null); do
  # Skip the line that defines project_root (we'll remove it)
  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD FIX: $(basename "$rscript")"
  else
    # Replace usage
    sed -i 's/project_root/ROOT_DIR/g' "$rscript"
    # Remove the now-redundant assignment (paths.R provides ROOT_DIR)
    sed -i '/^ROOT_DIR <- normalizePath(file\.path(script_dir/d' "$rscript"
    CHANGED=$((CHANGED + 1))
    echo "  Fixed: $(basename "$rscript")"
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
# 3. Replace output_dir = script_dir/outputs → RESULTS_DIR/intermediate/section
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- Step 3: Redirect output_dir from code/ to results/ ---"

for rscript in $(grep -rl 'file\.path(script_dir, "outputs")' "$CODE" --include="*.R" 2>/dev/null); do
  section=$(basename "$(dirname "$rscript")")
  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD FIX: $(basename "$rscript") → RESULTS_DIR/intermediate/${section}_outputs"
  else
    sed -i "s|file\.path(script_dir, \"outputs\")|file.path(RESULTS_DIR, \"intermediate\", \"${section}_outputs\")|g" "$rscript"
    CHANGED=$((CHANGED + 1))
    echo "  Fixed: $(basename "$rscript") → ${section}_outputs"
  fi
done

# Also fix figures_dir patterns
for rscript in $(grep -rl 'file\.path(script_dir, "figures")' "$CODE" --include="*.R" 2>/dev/null); do
  section=$(basename "$(dirname "$rscript")")
  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD FIX figures_dir: $(basename "$rscript")"
  else
    sed -i "s|file\.path(script_dir, \"figures\")|file.path(RESULTS_DIR, \"figures\", \"${section}\")|g" "$rscript"
    CHANGED=$((CHANGED + 1))
    echo "  Fixed figures_dir: $(basename "$rscript")"
  fi
done

# Fix fig_dir patterns that use project_root/figures/by_analysis/
for rscript in $(grep -rl 'figures/by_analysis' "$CODE" --include="*.R" 2>/dev/null); do
  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD FIX fig_dir: $(basename "$rscript")"
  else
    sed -i 's|file\.path(ROOT_DIR, "figures/by_analysis/\([^"]*\)")|file.path(RESULTS_DIR, "figures", "\1")|g' "$rscript"
    CHANGED=$((CHANGED + 1))
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
# 4. Replace data_dir patterns that use project_root/data/
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- Step 4: Replace data_dir → DATA_DIR ---"

# project_root/data/human_bulk_rnaseq → DATA_DIR/human_bulk_rnaseq
for rscript in $(grep -rl 'ROOT_DIR, "data/' "$CODE" --include="*.R" 2>/dev/null); do
  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD FIX data_dir: $(basename "$rscript")"
  else
    sed -i 's|file\.path(ROOT_DIR, "data/\([^"]*\)")|file.path(DATA_DIR, "\1")|g' "$rscript"
    CHANGED=$((CHANGED + 1))
    echo "  Fixed data_dir: $(basename "$rscript")"
  fi
done

# results/corrected/ → RESULTS_DIR/corrected/
for rscript in $(grep -rl 'ROOT_DIR, "results/' "$CODE" --include="*.R" 2>/dev/null); do
  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD FIX results_dir: $(basename "$rscript")"
  else
    sed -i 's|file\.path(ROOT_DIR, "results/\([^"]*\)")|file.path(RESULTS_DIR, "\1")|g' "$rscript"
    CHANGED=$((CHANGED + 1))
    echo "  Fixed results_dir: $(basename "$rscript")"
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
# 5. Arial font paths → system font fallback
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- Step 5: Font path fallbacks ---"

ARIAL_HPC='/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Arial.ttf'
LIBERATION='/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'

# R: font_add("Arial", "/ix1/.../Arial.ttf") → try/fallback
for rscript in $(grep -rl 'font_add.*Arial.*ix1' "$CODE" --include="*.R" 2>/dev/null); do
  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD FIX font: $(basename "$rscript")"
  else
    sed -i "s|font_add(\"Arial\", \"$ARIAL_HPC\")|tryCatch(font_add(\"Arial\", \"$LIBERATION\"), error = function(e) message(\"Arial font not available, using default\"))|g" "$rscript"
    CHANGED=$((CHANGED + 1))
    echo "  Fixed font: $(basename "$rscript")"
  fi
done

# Python: Replace HPC Arial paths AND clean up any stale fallback blocks
for pyscript in $(grep -rl "ix1.*Arial\|ARIAL.*ix1" "$CODE" --include="*.py" 2>/dev/null); do
  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD FIX font: $(basename "$pyscript")"
  else
    # Replace the hardcoded path assignment with a clean fallback
    sed -i "s|ARIAL_PATH = '$ARIAL_HPC'|ARIAL_PATH = next((f for f in ['$LIBERATION', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'] if os.path.exists(f)), None)|g" "$pyscript"
    sed -i "s|ARIAL_FONT_PATH = '$ARIAL_HPC'|ARIAL_FONT_PATH = next((f for f in ['$LIBERATION', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'] if os.path.exists(f)), None)|g" "$pyscript"
    sed -i "s|arial_path = '$ARIAL_HPC'|arial_path = next((f for f in ['$LIBERATION', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'] if os.path.exists(f)), None)|g" "$pyscript"
    # Remove ALL stale fallback lines (indented or not) that would break when value is already None
    sed -i '/^[[:space:]]*if not os\.path\.exists(ARIAL/d' "$pyscript"
    sed -i '/^[[:space:]]*if not os\.path\.exists(arial_path)/d' "$pyscript"
    sed -i '/Fallback to system Liberation/d' "$pyscript"
    # Remove orphaned indented reassignment lines left by previous fallback blocks
    sed -i '/^[[:space:]]*ARIAL_FONT_PATH = .*LiberationSans/d' "$pyscript"
    sed -i '/^[[:space:]]*ARIAL_PATH = .*LiberationSans/d' "$pyscript"
    sed -i '/^[[:space:]]*ARIAL_FONT_PATH = None/d' "$pyscript"
    sed -i '/^[[:space:]]*ARIAL_PATH = None.*matplotlib/d' "$pyscript"
    CHANGED=$((CHANGED + 1))
    echo "  Fixed font: $(basename "$pyscript")"
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
# 6. Section-specific hardcoded paths
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- Step 6: Section-specific path fixes ---"

# 01: Remove raw/ and external/ subdirectory prefixes (data asset has files at top level)
for f in $(find "$CODE/01_human_bulk_rnaseq" -name "*.R" 2>/dev/null); do
  if grep -q '"raw/' "$f" 2>/dev/null || grep -q '"external/' "$f" 2>/dev/null; then
    if [ "$DRY_RUN" = "1" ]; then
      echo "  WOULD strip raw/external subdirs: $(basename "$f")"
    else
      sed -i 's|, "raw/|, "|g' "$f"
      sed -i 's|, "external/|, "|g' "$f"
      CHANGED=$((CHANGED + 1))
      echo "  Stripped raw/external subdirs: $(basename "$f")"
    fi
  fi
done

# 02: Remove raw/ subdirectory prefix for rat_snrnaseq data paths
for f in $(find "$CODE/02_rat_snrnaseq" -name "*.R" 2>/dev/null); do
  if grep -q 'rat_snrnaseq/raw' "$f" 2>/dev/null || grep -q '"raw/' "$f" 2>/dev/null; then
    if [ "$DRY_RUN" = "1" ]; then
      echo "  WOULD strip raw subdir: $(basename "$f")"
    else
      sed -i 's|"rat_snrnaseq/raw"|"rat_snrnaseq"|g' "$f"
      sed -i 's|, "raw/|, "|g' "$f"
      CHANGED=$((CHANGED + 1))
      echo "  Stripped raw subdir: $(basename "$f")"
    fi
  fi
done

# 07: Fix _config.py SAMPLE_PATHS to use flat data structure (OS01/filtered_feature_bc_matrix/)
if [ -f "$CODE/07_organoid_single_cell/_config.py" ]; then
  # Replace pool1_flex/outs/per_sample_outs/OS01/count → OS01
  if grep -q 'pool1_flex/outs/per_sample_outs' "$CODE/07_organoid_single_cell/_config.py" 2>/dev/null; then
    if [ "$DRY_RUN" = "1" ]; then
      echo "  WOULD fix _config.py SAMPLE_PATHS to flat structure"
    else
      sed -i 's|pool1_cellranger.*per_sample_outs|_NEIL_PROJECT|g; s|pool2_cellranger.*per_sample_outs|_NEIL_PROJECT|g' "$CODE/07_organoid_single_cell/_config.py"
      CHANGED=$((CHANGED + 1))
      echo "  Fixed _config.py SAMPLE_PATHS to flat structure"
    fi
  fi
  # Fix sample path construction in _config.py to use OS01/filtered_feature_bc_matrix
  if grep -q 'sample_filtered_feature_bc_matrix' "$CODE/07_organoid_single_cell/_config.py" 2>/dev/null; then
    sed -i 's|/count/sample_filtered_feature_bc_matrix|/filtered_feature_bc_matrix|g' "$CODE/07_organoid_single_cell/_config.py"
    echo "  Fixed _config.py sample matrix path format"
  fi
fi

# 06: Ensure spatial scripts use DATA_DIR not relative 'data/' path
for f in "$CODE/06_spatial_biopsies/01_immune_secretion.py" "$CODE/06_spatial_biopsies/02_immune_pathways.py"; do
  if [ -f "$f" ] && grep -q "open('data/" "$f" 2>/dev/null; then
    if [ "$DRY_RUN" = "1" ]; then
      echo "  WOULD fix relative data/ path: $(basename "$f")"
    else
      # These should already use DATA_DIR from paths.py — if not, the paths import injection handles it
      echo "  NOTE: $(basename "$f") still uses relative data/ paths — needs manual fix if paths.py not imported"
    fi
  fi
done

# 01/mica/07: MICA internal functions
MICA_HPC='/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/MICA/R/00_internal_functions.R'
do_sed "$CODE/01_human_bulk_rnaseq/mica/07_statistical_audit.R" \
  "$MICA_HPC" \
  'file.path(CODE_DIR, "lib", "MICA", "R", "00_internal_functions.R")'

# 02: seurat obj hardcoded path
SEURAT_HPC='/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen/analysis/02_rat_snrnaseq/outputs/seurat_annotated_sctype_cluster.rds'
for f in "$CODE/02_rat_snrnaseq/method_comparison.R" "$CODE/02_rat_snrnaseq/spot_check_markers.R"; do
  do_sed "$f" "$SEURAT_HPC" 'resolve_input("02_rat_snrnaseq", "seurat_annotated_sctype_cluster.rds")'
done

# 03_wes: VCF input path — replace entire line to avoid nested quote issues
if [ -f "$CODE/03_rat_wes/02_cosmic_signatures.py" ]; then
  sed -i 's|^VCF_INPUT = Path(".*")|VCF_INPUT = DATA_DIR / "rat_wes"|' "$CODE/03_rat_wes/02_cosmic_signatures.py"
  CHANGED=$((CHANGED + 1))
  echo "  Fixed: 02_cosmic_signatures.py VCF_INPUT"
fi

# 03_wes: 01_parse_vep.py has script-relative paths that need capsule paths
if [ -f "$CODE/03_rat_wes/01_parse_vep.py" ]; then
  if grep -q 'Path(__file__).parent / "outputs"' "$CODE/03_rat_wes/01_parse_vep.py" 2>/dev/null; then
    if [ "$DRY_RUN" = "1" ]; then
      echo "  WOULD FIX: 01_parse_vep.py VEP_DIR and OUTPUT_DIR"
    else
      sed -i 's|VEP_DIR = Path(__file__).parent / "outputs" / "vep"|VEP_DIR = DATA_DIR / "rat_wes" / "vep"|' "$CODE/03_rat_wes/01_parse_vep.py"
      sed -i 's|OUTPUT_DIR = Path(__file__).parent / "outputs"|OUTPUT_DIR = RESULTS_DIR / "intermediate" / "03_rat_wes_outputs"|' "$CODE/03_rat_wes/01_parse_vep.py"
      # Add paths import if not present
      if ! grep -q 'from paths import' "$CODE/03_rat_wes/01_parse_vep.py" 2>/dev/null; then
        sed -i '/^from pathlib import Path/a\import sys; sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "config"))\nfrom paths import DATA_DIR, RESULTS_DIR' "$CODE/03_rat_wes/01_parse_vep.py"
      fi
      CHANGED=$((CHANGED + 1))
      echo "  Fixed: 01_parse_vep.py VEP_DIR + OUTPUT_DIR"
    fi
  fi
fi

# 04: Xu atlas path — replace entire line to avoid nested quote issues
if [ -f "$CODE/04_human_scrnaseq/00b_preprocess_seurat.R" ]; then
  sed -i 's|^raw_data_dir <- ".*Xu_etal.*"|raw_data_dir <- file.path(DATA_DIR, "human_scrnaseq")|' "$CODE/04_human_scrnaseq/00b_preprocess_seurat.R"
  CHANGED=$((CHANGED + 1))
  echo "  Fixed: 00b_preprocess_seurat.R raw_data_dir"
fi

# 05: Published validation dir — replace entire line
if [ -f "$CODE/05_rat_bulk_rnaseq/05_validate_vs_published.R" ]; then
  sed -i 's|^published_dir <- Sys.getenv.*|published_dir <- file.path(DATA_DIR, "rat_bulk_rnaseq", "published_validation")|' "$CODE/05_rat_bulk_rnaseq/05_validate_vs_published.R"
  CHANGED=$((CHANGED + 1))
  echo "  Fixed: 05_validate_vs_published.R published_dir"
fi

# 07: Neil project path — replace entire line to avoid nested quote issues
if [ -f "$CODE/07_organoid_single_cell/_config.py" ]; then
  sed -i 's|^_NEIL_PROJECT = ".*"|_NEIL_PROJECT = DATA_DIR / "organoid_single_cell"|' "$CODE/07_organoid_single_cell/_config.py"
  CHANGED=$((CHANGED + 1))
  echo "  Fixed: _config.py _NEIL_PROJECT"
fi

# Comment out FASTQ paths in _config.py (Cell Ranger excluded from capsule)
if grep -q 'fastqs_pool1.*ix1' "$CODE/07_organoid_single_cell/_config.py" 2>/dev/null; then
  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD comment out FASTQ paths in _config.py"
  else
    sed -i 's|"fastqs_pool1": "/ix1/|# "fastqs_pool1": "/ix1/|' "$CODE/07_organoid_single_cell/_config.py"
    sed -i 's|"fastqs_pool2": "/ix1/|# "fastqs_pool2": "/ix1/|' "$CODE/07_organoid_single_cell/_config.py"
    CHANGED=$((CHANGED + 1))
    echo "  Commented out FASTQ paths in _config.py"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# 7. Python scripts: inject paths.py import where needed
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- Step 7: Inject paths.py imports into Python scripts ---"

for pyscript in $(grep -rl 'DATA_DIR\|RESULTS_DIR\|resolve_input' "$CODE" --include="*.py" 2>/dev/null); do
  # Skip config/ and already-imported files
  if echo "$pyscript" | grep -q "/config/"; then continue; fi
  if grep -q "from paths import" "$pyscript" 2>/dev/null; then continue; fi

  if [ "$DRY_RUN" = "1" ]; then
    echo "  WOULD INJECT paths.py import: $(basename "$pyscript")"
  else
    # Insert after the last import line
    sed -i '0,/^import \|^from /!b; /^import \|^from /{h;$!{n;/^import \|^from /!{x;a\
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "config"))\
from paths import DATA_DIR, RESULTS_DIR, resolve_input, save_intermediate, save_figure
x}};}' "$pyscript" 2>/dev/null || true
    CHANGED=$((CHANGED + 1))
    echo "  Injected paths.py: $(basename "$pyscript")"
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
# 8. Verify no HPC paths remain in active code
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "--- Step 8: Verification ---"

REMAINING=$(grep -rn '/ix1/\|/ihome/\|/bgfs/' "$CODE" --include="*.R" --include="*.py" | grep -v '^.*:#' | wc -l)
COMMENTS=$(grep -rn '/ix1/\|/ihome/\|/bgfs/' "$CODE" --include="*.R" --include="*.py" | grep '^.*:#' | wc -l)

echo "  HPC paths in active code: $REMAINING"
echo "  HPC paths in comments: $COMMENTS (OK)"

if [ "$REMAINING" -gt 0 ]; then
  echo ""
  echo "  WARNING: Remaining HPC paths in active code:"
  grep -rn '/ix1/\|/ihome/\|/bgfs/' "$CODE" --include="*.R" --include="*.py" | grep -v '^.*:#'
fi

echo ""
echo "================================================================"
echo "  Files modified: $CHANGED"
if [ "$DRY_RUN" = "1" ]; then
  echo "  (dry-run — nothing was written)"
fi
echo "================================================================"
