#!/usr/bin/env python
"""Bundle all analysis figures with a manifest mapping to PPT slides.

Creates a flat export directory with:
  - All figure PNGs/SVGs copied with descriptive names
  - manifest.csv mapping each figure to its analysis section and PPT slide

Usage:
    python scripts/bundle_figures.py [--output-dir figures/bundle]
"""
import shutil
import csv
import sys
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = PROJECT_ROOT / "figures" / "bundle"

# Figure manifest: (source_path_relative, tag, ppt_slide_note)
# tag = which manuscript/PPT figure this corresponds to
FIGURE_MANIFEST = [
    # === Section 03: Rat WES ===
    ("analysis/03_rat_wes/figures/cosmic_signatures.png",
     "Fig2B", "PPT Slide 4 — COSMIC Signatures (PASS)"),
    ("analysis/03_rat_wes/figures/cosmic_signatures.svg",
     "Fig2B", "PPT Slide 4 — COSMIC Signatures (SVG)"),
    ("analysis/03_rat_wes/figures/oncoplot.png",
     "Fig2_supp", "Oncoplot — supporting"),
    ("analysis/03_rat_wes/figures/oncoplot.svg",
     "Fig2_supp", "Oncoplot — supporting (SVG)"),

    # === Section 04: Human scRNA-seq ===
    ("analysis/04_human_scrnaseq/figures/gsva_heatmap.png",
     "Fig4E", "PPT Slide 6 — Estrogen Pathway Heatmap (FIXED)"),
    ("analysis/04_human_scrnaseq/figures/gsva_heatmap.svg",
     "Fig4E", "Estrogen Pathway Heatmap (SVG)"),
    ("analysis/04_human_scrnaseq/figures/fig7bc_pathway_heatmaps.png",
     "Fig7BC", "PPT Slides 8-9 — HALLMARK/BIOCARTA Pathway Heatmaps (FIXED)"),
    ("analysis/04_human_scrnaseq/figures/fig7bc_pathway_heatmaps.svg",
     "Fig7BC", "HALLMARK/BIOCARTA Pathway Heatmaps (SVG)"),
    ("analysis/04_human_scrnaseq/figures/fig7d_cellphonedb_dotplot.png",
     "Fig7D", "PPT Slide 10 — CellPhoneDB Communication (FIXED)"),
    ("analysis/04_human_scrnaseq/figures/fig7d_cellphonedb_dotplot.svg",
     "Fig7D", "CellPhoneDB Communication (SVG)"),
    ("analysis/04_human_scrnaseq/figures/fig7d_cellphonedb_dotplot_curated.png",
     "Fig7D_curated", "Curated L-R pairs subset"),
    ("analysis/04_human_scrnaseq/figures/fig7d_cellchat_dotplot.png",
     "Fig7D_alt", "Alternative CellChat version"),
    ("analysis/04_human_scrnaseq/figures/fig7c_pathway_heatmap.png",
     "Fig7C_macro", "Macrophage-specific pathway heatmap"),
    ("analysis/04_human_scrnaseq/figures/fig7c_pathway_barplot.png",
     "Fig7C_bar", "Pathway enrichment barplot"),
    ("analysis/04_human_scrnaseq/figures/fig7b_macrophage_volcano.png",
     "Fig7B_volcano", "Macrophage DEG volcano"),
    ("analysis/04_human_scrnaseq/figures/fig7b_macrophage_heatmap.png",
     "Fig7B_heatmap", "Macrophage DEG heatmap"),
    ("analysis/04_human_scrnaseq/figures/fig7b_celltype_deg_barplot.png",
     "Fig7B_bar", "Cell type DEG barplot"),
    ("analysis/04_human_scrnaseq/figures/umap_celltypes.png",
     "Fig_supp_umap", "Cell type UMAP"),
    ("analysis/04_human_scrnaseq/figures/progeny_heatmap.png",
     "Fig_supp_progeny", "PROGENy pathway heatmap"),
    ("analysis/04_human_scrnaseq/figures/feature_markers.png",
     "Fig_supp_markers", "Feature marker plots"),

    # === Section 06: Spatial biopsies ===
    ("analysis/06_spatial_biopsies/figures/immune_pathways/dotplots/prerank_summary_dotplot.png",
     "Fig_spatial_summary", "Immune pathway summary dotplot"),
    ("analysis/06_spatial_biopsies/figures/immune_pathways/dotplots/prerank_summary_dotplot.svg",
     "Fig_spatial_summary", "Immune pathway summary dotplot (SVG)"),

    # === Section 07: Organoid single cell (ICI→fulv updated) ===
    # These will be regenerated — include paths for when they exist
    ("analysis/07_organoid_single_cell/figures/single_cell_pathways/heatmap_effect_sizes.png",
     "Fig_organoid_pathways", "Treatment pathway effect sizes"),
    ("analysis/07_organoid_single_cell/figures/single_cell_pathways/hypothesis_summary.png",
     "Fig_organoid_hypothesis", "Hypothesis test summary"),
    ("analysis/07_organoid_single_cell/figures/inhibitor_mechanism/signature_projection.png",
     "Fig_organoid_signatures", "E1/E2 signature projection"),
    ("analysis/07_organoid_single_cell/figures/inhibitor_mechanism/fulv_calibration_kde.png",
     "Fig_organoid_fulv_cal", "Fulvestrant calibration KDE"),
    ("analysis/07_organoid_single_cell/figures/inhibitor_mechanism/subprogram_effect_heatmap.png",
     "Fig_organoid_subprog", "ER sub-program effect heatmap"),
    ("analysis/07_organoid_single_cell/figures/inhibitor_mechanism/dose_response_hsd17b7.png",
     "Fig_organoid_dose", "HSD17B7 dose-response"),
    ("analysis/07_organoid_single_cell/figures/proliferation/panel1_score_violins.png",
     "Fig_organoid_prolif1", "Proliferation score violins"),
    ("analysis/07_organoid_single_cell/figures/proliferation/panel2_cycling_fraction.png",
     "Fig_organoid_prolif2", "Cycling fraction barplot"),
    ("analysis/07_organoid_single_cell/figures/proliferation/panel6_cluster_cycling_heatmap.png",
     "Fig_organoid_prolif6", "Per-cluster cycling heatmap"),
    ("analysis/07_organoid_single_cell/figures/proliferation/panel7_prolif_er_coupling.png",
     "Fig_organoid_prolif7", "Proliferation-ER coupling"),
    ("analysis/07_organoid_single_cell/figures/heterogeneity/01_umap_density_per_condition.png",
     "Fig_organoid_het1", "UMAP density per condition"),
    ("analysis/07_organoid_single_cell/figures/heterogeneity/05_biological_axis.png",
     "Fig_organoid_het5", "Biological axis PCA"),
    ("analysis/07_organoid_single_cell/figures/cell_cycle/phase_stacked_bar.png",
     "Fig_organoid_cc", "Cell cycle phase stacked bar"),
]


def main():
    output_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUTPUT
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{datetime.now():%Y-%m-%d %H:%M}] Bundling figures to {output_dir}")

    manifest_rows = []
    copied = 0
    missing = 0

    for rel_path, tag, note in FIGURE_MANIFEST:
        src = PROJECT_ROOT / rel_path
        if src.exists():
            # Flat name: tag + original extension
            suffix = src.suffix
            dest_name = f"{tag}_{src.stem}{suffix}"
            dest = output_dir / dest_name
            shutil.copy2(src, dest)
            manifest_rows.append({
                "bundle_filename": dest_name,
                "tag": tag,
                "source": rel_path,
                "ppt_note": note,
                "size_kb": round(src.stat().st_size / 1024, 1),
            })
            copied += 1
        else:
            manifest_rows.append({
                "bundle_filename": "(MISSING)",
                "tag": tag,
                "source": rel_path,
                "ppt_note": note + " [NOT YET GENERATED]",
                "size_kb": 0,
            })
            missing += 1

    # Write manifest CSV
    manifest_path = output_dir / "manifest.csv"
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["bundle_filename", "tag", "source", "ppt_note", "size_kb"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"  Copied: {copied} figures")
    print(f"  Missing (not yet regenerated): {missing}")
    print(f"  Manifest: {manifest_path}")
    print(f"  Done.")


if __name__ == "__main__":
    main()
