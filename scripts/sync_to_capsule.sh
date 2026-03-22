#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────
# sync_to_capsule.sh
#
# One-directional sync: ERpBRCA_OlderWomen/analysis/ → ERpBRCA_CodeOcean/code/
# The main repo is the source of truth. Never edit capsule code directly.
# ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SOURCE_DIR="${REPO_ROOT}/analysis"
DEFAULT_TARGET="$(cd "${REPO_ROOT}/.." && pwd)/ERpBRCA_CodeOcean/code"

# ── Defaults ────────────────────────────────────────────────────────
TARGET_DIR="${DEFAULT_TARGET}"
DRY_RUN=0
DIFF_MODE=0
SECTION_FILTER=""

# ── Per-section inclusion lists ─────────────────────────────────────
declare -A INCLUDE

INCLUDE[01_human_bulk_rnaseq]="01_preprocess.R 02_run_gsva.R 03_run_progeny.R 04_correlations.R 05_visualize.R"
INCLUDE[01_human_bulk_rnaseq/mica]="01_prep_data.R 02_run_gsva.R 03_run_mica.R 04_fig4_estrogen.R 05_fig5_immune.R 07_statistical_audit.R"
INCLUDE[02_rat_snrnaseq]="01_qc_filter.R 02_normalize_integrate.R 03_cluster_annotate.R 03b_sctype_annotate.R 03c_sctype_cluster.R 08_differential_expression.R 09_differential_abundance.R method_comparison.R spot_check_markers.R"
INCLUDE[03_comparison]="01_generate_comparison_report.R"
INCLUDE[03_rat_wes]="01_parse_vep.py 02_cosmic_signatures.py 03_generate_oncoplot.py"
INCLUDE[04_human_scrnaseq]="00b_preprocess_seurat.R 01_load_subset_data.R 03_cell_type_annotation.R 04_cell_fractions.R 05_gene_expression_violin.R 06_run_gsva.R 07_run_progeny.R 08_run_wcsea.R 09_cellphonedb_prep.R 10_load_macrophage_seurat.R 11_macrophage_deg_analysis.R 12_macrophage_pathway_enrichment.R 12_enrichr_pathway_analysis.py 13_macrophage_cellphonedb.R 14_multicelltype_deg.R 15_multicelltype_pathway.R 16_cellchat_analysis.R 16_run_cellphonedb.py 17_cellphonedb_dotplot.R"
INCLUDE[05_rat_bulk_rnaseq]="03_deseq2.R 04_pam50_subtyping.R 05_validate_vs_rahul.R"
INCLUDE[06_spatial_biopsies]="01_immune_secretion.py 02_immune_pathways.py figure_config.py utils.py"
INCLUDE[07_organoid_single_cell]="02_qc.py 03_preprocess.py 04_pseudobulk_de.py 05_pathways.py 06_cell_cycle.py 07_single_cell_pathways.py 08_inhibitor_mechanism_exploration.py 09_pathway_enrichment_analysis.py 10_heterogeneity_analysis.py 11_proliferation_analysis.py _config.py _figure_config.py"

# ── Files that must never be copied ─────────────────────────────────
EXCLUDE_PATTERNS="run.sh __pycache__ .pipeline_markers sbatch_ submit_ run_analysis.sh run_analysis.sbatch"

# ── Usage ───────────────────────────────────────────────────────────
usage() {
    cat <<'USAGE'
Usage: sync_to_capsule.sh [OPTIONS]

Options:
  --section SECTION   Sync only the named section (e.g. 04_human_scrnaseq)
  --dry-run           Preview what would be copied without writing
  --diff              Show diff between source and existing capsule code
  --target PATH       Override capsule code directory
  -h, --help          Show this help message

Examples:
  ./scripts/sync_to_capsule.sh
  ./scripts/sync_to_capsule.sh --section 04_human_scrnaseq
  ./scripts/sync_to_capsule.sh --dry-run
  ./scripts/sync_to_capsule.sh --diff
  ./scripts/sync_to_capsule.sh --target /path/to/capsule/code
USAGE
    exit 0
}

# ── Parse arguments ─────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --section)   SECTION_FILTER="$2"; shift 2 ;;
        --dry-run)   DRY_RUN=1; shift ;;
        --diff)      DIFF_MODE=1; shift ;;
        --target)    TARGET_DIR="$2"; shift 2 ;;
        -h|--help)   usage ;;
        *)           echo "ERROR: Unknown option: $1" >&2; usage ;;
    esac
done

# ── Validate paths ──────────────────────────────────────────────────
if [[ ! -d "${SOURCE_DIR}" ]]; then
    echo "ERROR: Source directory not found: ${SOURCE_DIR}" >&2
    exit 1
fi

if [[ "${DIFF_MODE}" -eq 0 && "${DRY_RUN}" -eq 0 && ! -d "${TARGET_DIR}" ]]; then
    echo "ERROR: Target directory not found: ${TARGET_DIR}" >&2
    echo "       Use --target to specify a different capsule code path." >&2
    exit 1
fi

# ── Helpers ─────────────────────────────────────────────────────────
is_excluded() {
    local filename="$1"
    for pat in ${EXCLUDE_PATTERNS}; do
        case "${filename}" in
            ${pat}*) return 0 ;;
        esac
    done
    return 1
}

git_file_info() {
    local filepath="$1"
    # Returns "short_sha date" or "unknown unknown" if not tracked
    local info
    info="$(git -C "${REPO_ROOT}" log -1 --format="%h %ai" -- "${filepath}" 2>/dev/null || true)"
    if [[ -z "${info}" ]]; then
        echo "untracked $(date -r "${filepath}" '+%Y-%m-%d %H:%M:%S %z' 2>/dev/null || echo 'unknown')"
    else
        echo "${info}"
    fi
}

# ── Counters ────────────────────────────────────────────────────────
COPIED=0
SKIPPED=0
ERRORS=0
MANIFEST_LINES=()

# ── Determine which sections to process ─────────────────────────────
SECTIONS_TO_PROCESS=()
for section in "${!INCLUDE[@]}"; do
    if [[ -n "${SECTION_FILTER}" ]]; then
        # Match exact section or any sub-path (e.g. --section 01_human_bulk_rnaseq matches 01_human_bulk_rnaseq/mica)
        if [[ "${section}" == "${SECTION_FILTER}" || "${section}" == "${SECTION_FILTER}/"* || "${SECTION_FILTER}" == "${section}/"* ]]; then
            SECTIONS_TO_PROCESS+=("${section}")
        fi
    else
        SECTIONS_TO_PROCESS+=("${section}")
    fi
done

if [[ ${#SECTIONS_TO_PROCESS[@]} -eq 0 ]]; then
    if [[ -n "${SECTION_FILTER}" ]]; then
        echo "ERROR: No matching section for filter '${SECTION_FILTER}'" >&2
        echo "Available sections:" >&2
        for s in "${!INCLUDE[@]}"; do echo "  ${s}" >&2; done
        exit 1
    fi
    echo "ERROR: No sections defined." >&2
    exit 1
fi

# Sort sections for deterministic output
IFS=$'\n' SECTIONS_TO_PROCESS=($(sort <<<"${SECTIONS_TO_PROCESS[*]}")); unset IFS

# ── Main sync loop ──────────────────────────────────────────────────
echo "================================================================"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "  DRY RUN — no files will be written"
elif [[ "${DIFF_MODE}" -eq 1 ]]; then
    echo "  DIFF MODE — showing differences"
else
    echo "  SYNC: analysis/ -> $(basename "$(dirname "${TARGET_DIR}")")/code/"
fi
echo "  Source:  ${SOURCE_DIR}"
echo "  Target:  ${TARGET_DIR}"
echo "================================================================"
echo ""

for section in "${SECTIONS_TO_PROCESS[@]}"; do
    file_list="${INCLUDE[${section}]}"
    echo "── ${section} ──"

    for file in ${file_list}; do
        src="${SOURCE_DIR}/${section}/${file}"
        dst="${TARGET_DIR}/${section}/${file}"
        rel_src="analysis/${section}/${file}"
        rel_dst="code/${section}/${file}"

        # Safety: skip excluded files (should not appear in INCLUDE, but guard)
        if is_excluded "$(basename "${file}")"; then
            echo "  SKIP (excluded)  ${rel_src}"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        # Check source exists
        if [[ ! -f "${src}" ]]; then
            echo "  MISSING SOURCE   ${rel_src}"
            ERRORS=$((ERRORS + 1))
            continue
        fi

        # ── Diff mode ───────────────────────────────────────────────
        if [[ "${DIFF_MODE}" -eq 1 ]]; then
            if [[ -f "${dst}" ]]; then
                if diff -q "${src}" "${dst}" >/dev/null 2>&1; then
                    echo "  IDENTICAL  ${rel_dst}"
                else
                    echo "  DIFFERS    ${rel_dst}"
                    diff --color=auto -u "${dst}" "${src}" || true
                    echo ""
                fi
            else
                echo "  NEW        ${rel_dst}  (not yet in capsule)"
            fi
            continue
        fi

        # ── Dry-run mode ────────────────────────────────────────────
        if [[ "${DRY_RUN}" -eq 1 ]]; then
            echo "  COPY  ${rel_src}  ->  ${rel_dst}"
            COPIED=$((COPIED + 1))
            # Still collect manifest info
            info="$(git_file_info "${rel_src}")"
            sha="$(echo "${info}" | awk '{print $1}')"
            mod_date="$(echo "${info}" | awk '{print $2, $3}')"
            MANIFEST_LINES+=("| ${rel_dst} | ${rel_src} | ${sha} | ${mod_date} |")
            continue
        fi

        # ── Actual copy ─────────────────────────────────────────────
        mkdir -p "$(dirname "${dst}")"
        cp "${src}" "${dst}"
        echo "  COPY  ${rel_src}  ->  ${rel_dst}"
        COPIED=$((COPIED + 1))

        # Collect manifest info
        info="$(git_file_info "${rel_src}")"
        sha="$(echo "${info}" | awk '{print $1}')"
        mod_date="$(echo "${info}" | awk '{print $2, $3}')"
        MANIFEST_LINES+=("| ${rel_dst} | ${rel_src} | ${sha} | ${mod_date} |")
    done
    echo ""
done

# ── Generate MANIFEST.md ────────────────────────────────────────────
if [[ "${DIFF_MODE}" -eq 0 && "${DRY_RUN}" -eq 0 && ${#MANIFEST_LINES[@]} -gt 0 ]]; then
    MANIFEST="${TARGET_DIR}/MANIFEST.md"
    {
        echo "# CodeOcean Capsule Manifest"
        echo ""
        echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')"
        echo ""
        echo "| Capsule Path | Source Path | Git SHA | Modified Date |"
        echo "|---|---|---|---|"
        for line in "${MANIFEST_LINES[@]}"; do
            echo "${line}"
        done
    } > "${MANIFEST}"
    echo "── Generated ${MANIFEST} ──"
    echo ""
fi

# ── Validation: check run.sh references ─────────────────────────────
if [[ "${DIFF_MODE}" -eq 0 ]]; then
    echo "── Validation ──"
    VALIDATION_WARNINGS=0

    # Get the top-level section names (without /mica etc.)
    declare -A TOP_SECTIONS
    for section in "${SECTIONS_TO_PROCESS[@]}"; do
        top="${section%%/*}"
        TOP_SECTIONS["${top}"]=1
    done

    for top in $(echo "${!TOP_SECTIONS[@]}" | tr ' ' '\n' | sort); do
        run_sh="${TARGET_DIR}/${top}/run.sh"
        if [[ ! -f "${run_sh}" ]]; then
            continue
        fi

        # Extract script filenames referenced in run.sh (Rscript/python calls)
        while IFS= read -r raw_ref; do
            # Skip empty lines
            [[ -z "${raw_ref}" ]] && continue
            # Strip quotes, $SECTION_DIR/ prefix, and similar variable expansions
            referenced="${raw_ref//\"/}"
            referenced="${referenced//\'/}"
            referenced="${referenced#\$SECTION_DIR/}"
            referenced="${referenced#\$\{SECTION_DIR\}/}"
            referenced="${referenced#./}"
            [[ -z "${referenced}" ]] && continue
            # Only check actual script files (skip variable-only references)
            [[ "${referenced}" == \$* ]] && continue
            # Check if the referenced file exists in capsule
            if [[ ! -f "${TARGET_DIR}/${top}/${referenced}" ]]; then
                echo "  WARN: ${top}/run.sh references '${referenced}' but it is missing from capsule"
                VALIDATION_WARNINGS=$((VALIDATION_WARNINGS + 1))
            fi
        done < <(grep -oP '(?:Rscript|python3?)\s+\K[^\s;|&]+' "${run_sh}" 2>/dev/null || true)
    done

    if [[ "${VALIDATION_WARNINGS}" -eq 0 ]]; then
        echo "  All run.sh references satisfied."
    fi
    echo ""
fi

# ── Summary ─────────────────────────────────────────────────────────
echo "================================================================"
if [[ "${DIFF_MODE}" -eq 1 ]]; then
    echo "  Diff complete."
else
    echo "  Files copied:  ${COPIED}"
    echo "  Files skipped: ${SKIPPED}"
    echo "  Missing source: ${ERRORS}"
    echo "  Total:         $((COPIED + SKIPPED + ERRORS))"
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        echo "  (dry-run — nothing was written)"
    fi
fi
echo "================================================================"
