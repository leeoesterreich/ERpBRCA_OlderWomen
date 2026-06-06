#!/usr/bin/env bash
# Build EPS bundle for Neil - accepted Nature Aging manuscript figures
# Outputs: bundle/eps_for_neil/*.eps + manifest.csv + README.md + zip
set -euo pipefail

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen"
OD="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Onedrive_data/A.V. Lee Lab's files - Alex Chang/1. Projects/Neil Aging Project"
INK="/software/rhel9/manual/install/inkscape/1.4.2/Inkscape-ebf0e94-x86_64.AppImage"
OUT="$REPO/bundle/eps_for_neil"
DATE=$(date -u +%Y-%m-%d)

# NOTE: Iteration-state driver. Some entries below are pre-final (REGENERATE flags
# refer to early stages; final bundle is produced by separate single-panel scripts
# in analysis/<section>/). Kept for archival reference.

mkdir -p "$OUT"
MANIFEST="$OUT/manifest.csv"
echo "panel,subpanel,source,output_eps,source_status,notes" > "$MANIFEST"

# Final-state entries (post-iteration)
# Each row: panel|subpanel|source-path|output|status|notes
ENTRIES=(
  "2A||/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/output_plots/oncoplot_BRCA.svg|2A.eps|READY|Rat WES oncoplot (NeilRatWES original)"
  "2E||$REPO/analysis/02_rat_snrnaseq/figures/fig2e_rat_pathway_heatmap.svg|2E.eps|READY|Rat HALLMARK GSVA per celltype x age (10_multicelltype_pathway.R)"
  "5D|full|$REPO/analysis/07_organoid_single_cell/figures/heterogeneity/04_cluster_composition.svg|5D_full.eps|READY|Composite (full 3-subpanel)"
  "5D|cropped|$REPO/analysis/07_organoid_single_cell/figures/heterogeneity/04_cluster_umap_only.svg|5D_cropped.eps|READY|UMAP only (matches caption)"
  "5E|full|$REPO/analysis/07_organoid_single_cell/figures/heterogeneity/06c_estrogen_overlap_umap.svg|5E_full.eps|READY|Composite (heatmap + UMAP)"
  "5E|cropped|$REPO/analysis/07_organoid_single_cell/figures/heterogeneity/06c_estrogen_overlap_heatmap_only.svg|5E_cropped.eps|READY|Heatmap only (matches caption)"
  "5F||$REPO/analysis/07_organoid_single_cell/figures/inhibitor_mechanism/subprogram_effect_heatmap.svg|5F.eps|READY|ER sub-program effect"
  "5G||$REPO/analysis/07_organoid_single_cell/figures/inhibitor_mechanism/fulv_calibration_kde.svg|5G.eps|READY|Pathway-score KDE"
  "5H|activated|$REPO/analysis/07_organoid_single_cell/figures/inhibitor_mechanism/enrichr_activated_by_hsd17b7i_(e1_context).svg|5H_activated.eps|READY|Activated bars"
  "5H|suppressed|$REPO/analysis/07_organoid_single_cell/figures/inhibitor_mechanism/enrichr_suppressed_by_hsd17b7i_(e1_context).svg|5H_suppressed.eps|READY|Suppressed bars"
  "5I||$REPO/analysis/07_organoid_single_cell/figures/inhibitor_mechanism/gsea_e1_vs_e1_hsd17b7i.svg|5I.eps|READY|GSEA dotplot"
  "7B||$REPO/analysis/04_human_scrnaseq/figures/fig7bc_pathway_heatmaps.svg|7B.eps|READY|Hallmark+BIOCARTA heatmaps"
  "7D||$REPO/analysis/04_human_scrnaseq/figures/fig7d_elderly_minus_young.svg|7D.eps|READY|Elderly-Young diff heatmap (15c_fig7d_elderly_minus_young.R)"
  "7G|P1|$REPO/analysis/06_spatial_biopsies/figures/fig7g/HCC22-088-P1-S1_immunosuppressive_macrophages_spatial.svg|7G_P1.eps|READY|CD163+ proportion spatial"
  "7G|P2|$REPO/analysis/06_spatial_biopsies/figures/fig7g/HCC22-088-P2-S1_immunosuppressive_macrophages_spatial.svg|7G_P2.eps|READY|CD163+ proportion spatial"
  "7G|P3|$REPO/analysis/06_spatial_biopsies/figures/fig7g/HCC22-088-P3-S1_A_immunosuppressive_macrophages_spatial.svg|7G_P3.eps|READY|CD163+ proportion spatial"
  "7G|P4|$REPO/analysis/06_spatial_biopsies/figures/fig7g/HCC22-088-P4-S1_immunosuppressive_macrophages_spatial.svg|7G_P4.eps|READY|CD163+ proportion spatial"
  "7G|P5|$REPO/analysis/06_spatial_biopsies/figures/fig7g/HCC22-088-P5-S1_immunosuppressive_macrophages_spatial.svg|7G_P5.eps|READY|CD163+ proportion spatial"
  "7G|P6|$REPO/analysis/06_spatial_biopsies/figures/fig7g/HCC22-088-P6-S1_immunosuppressive_macrophages_spatial.svg|7G_P6.eps|READY|CD163+ proportion spatial"
  "7H||$REPO/analysis/06_spatial_biopsies/figures/immune_pathways/dotplots/prerank_summary_dotplot.svg|7H.eps|READY|IL/TGF-beta pathway heatmap"
  "7I||$OD/neil_aging_biopsy_figures/Immunosuppressive_Macrophages_CD163plus/Immunosuppressive_Macrophages_CD163plus_sender_top_pathways.svg|7I_sender_top_pathways.eps|READY|COMMOT sender pathways"
  "8I||/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/figures/colocalization_analysis/mean_colocalization_heatmap_combined_macrophages.svg|8I.eps|READY|Cell-type co-localization"
)

convert_svg_to_eps() {
  local src="$1" dst="$2"
  "$INK" --appimage-extract-and-run \
    --export-type=eps --export-text-to-path=false \
    --export-filename="$dst" "$src" 2>&1 | grep -vE "GtkRecentManager|appimage_extracted|Failed to clean up cache" || true
}

for entry in "${ENTRIES[@]}"; do
  IFS='|' read -r panel sub src out status notes <<<"$entry"
  out_path="$OUT/$out"
  [ ! -f "$src" ] && { echo "MISSING SOURCE: $src" >&2; echo "$panel,$sub,$src,$out,MISSING,$notes" >> "$MANIFEST"; continue; }
  echo ">> $out"
  convert_svg_to_eps "$src" "$out_path"
  $INK --appimage-extract-and-run --export-type=pdf \
    --export-filename="$OUT/pdf_review/${out%.eps}.pdf" "$src" >/dev/null 2>&1
  echo "$panel,$sub,$src,$out,$status,$notes" >> "$MANIFEST"
done

ZIPNAME="$REPO/bundle/eps_for_neil_${DATE}.zip"
(cd "$REPO/bundle" && zip -qr "$(basename "$ZIPNAME")" "eps_for_neil")
echo "Bundle: $OUT"
echo "Zip:    $ZIPNAME"
