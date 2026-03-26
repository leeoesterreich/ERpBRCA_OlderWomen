#!/usr/bin/env python
# analysis/04_human_scrnaseq/12_enrichr_pathway_analysis.py
# ENRICHR pathway analysis via GSEApy on pseudobulk DEGs (Elderly vs Young)
#
# Runs ENRICHR on significant DEGs (FDR < 0.05) from each cell type
# using the multicelltype pseudobulk DEG results (script 14) and
# macrophage-specific DEG results (script 11).
#
# Inputs:
#   - outputs/multicelltype_deg_full/*.csv   (per-celltype DEG tables from DESeq2)
#   - outputs/macrophage_degs.csv            (macrophage-specific DEG table)
#
# Outputs:
#   - outputs/enrichr_results/               (per-celltype ENRICHR results)
#   - outputs/enrichr_summary.csv            (all cell types, top pathways)
#   - figures/fig7c_enrichr_heatmap.png       (heatmap: cell types x Hallmark pathways)
#   - figures/fig7c_enrichr_heatmap.pdf
#   - figures/fig7c_enrichr_macrophage.png    (macrophage-specific pathways)

import os
import sys
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

np.random.seed(12345)

# Directories
script_dir = Path(__file__).resolve().parent
output_dir = script_dir / "outputs"
figures_dir = script_dir / "figures"
deg_dir = output_dir / "multicelltype_deg_full"
enrichr_dir = output_dir / "enrichr_results"
enrichr_dir.mkdir(parents=True, exist_ok=True)
figures_dir.mkdir(parents=True, exist_ok=True)

# ENRICHR libraries to query (via API)
ENRICHR_LIBS = [
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "Reactome_Pathways_2024",
    "GO_Biological_Process_2023",
    "WikiPathway_2023_Human",
    "BioCarta_2016",
]

# Local GMT files (not on ENRICHR, run offline Fisher's exact via gseapy.enrich)
LOCAL_GMT_LIBS = {}

FDR_THRESHOLD = 0.05
# For macrophage-only (0 FDR<0.05 DEGs), use relaxed nominal p-value
NOMINAL_P_THRESHOLD = 0.01


def clean_gene_name(gene):
    """Fix tab-duplicated gene names (e.g., 'GEM\\tGEM' → 'GEM')."""
    if "\t" in str(gene):
        return str(gene).split("\t")[0]
    return str(gene)


def load_degs(filepath, fdr_thresh=FDR_THRESHOLD, use_nominal=False):
    """Load DEG table, clean gene names, return up/down gene lists."""
    df = pd.read_csv(filepath)
    df["gene"] = df["gene"].apply(clean_gene_name)

    if use_nominal:
        # Use nominal p-value for underpowered comparisons
        sig = df[df["pvalue"] < NOMINAL_P_THRESHOLD].copy()
        thresh_label = f"p<{NOMINAL_P_THRESHOLD}"
    else:
        sig = df[df["padj"] < fdr_thresh].copy()
        thresh_label = f"FDR<{fdr_thresh}"

    up = sig[sig["log2FoldChange"] > 0]["gene"].tolist()
    down = sig[sig["log2FoldChange"] < 0]["gene"].tolist()
    all_genes = df["gene"].tolist()

    return up, down, all_genes, thresh_label


def run_enrichr_on_genes(gene_list, description, outdir):
    """Run ENRICHR on a gene list, one library at a time for correct FDR.

    Running libraries individually ensures adjusted p-values are computed
    within each library rather than across all libraries combined.
    """
    if len(gene_list) == 0:
        return pd.DataFrame()

    all_results = []
    for lib in ENRICHR_LIBS:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=lib,
                organism="human",
                outdir=str(outdir),
                no_plot=True,
                cutoff=1.0,
            )
            res = enr.results.copy()
            res["description"] = description
            all_results.append(res)
        except Exception as e:
            print(f"    WARNING: ENRICHR failed for {description}/{lib}: {e}")

    # Run local GMT libraries (offline Fisher's exact test via gseapy.enrich)
    for lib_name, gmt_path in LOCAL_GMT_LIBS.items():
        if not gmt_path.exists():
            print(f"    WARNING: GMT file not found: {gmt_path}")
            continue
        try:
            enr = gp.enrich(
                gene_list=gene_list,
                gene_sets=str(gmt_path),
                background=20000,  # approx human protein-coding genes
                outdir=str(outdir),
                no_plot=True,
                cutoff=1.0,
            )
            res = enr.results.copy()
            res["Gene_set"] = lib_name
            res["description"] = description
            all_results.append(res)
        except Exception as e:
            print(f"    WARNING: enrich() failed for {description}/{lib_name}: {e}")

    if len(all_results) == 0:
        return pd.DataFrame()
    return pd.concat(all_results, ignore_index=True)


def main():
    print("=== ENRICHR Pathway Analysis via GSEApy ===\n")

    # -------------------------------------------------------------------------
    # Step 1: Multi-celltype DEGs (from script 14)
    # -------------------------------------------------------------------------
    print("Step 1: Loading multicelltype DEG results...")

    all_enrichr = []
    deg_files = sorted(deg_dir.glob("*_degs.csv"))
    print(f"  Found {len(deg_files)} cell type DEG files\n")

    for deg_file in deg_files:
        ct_name = deg_file.stem.replace("_degs", "")
        up, down, all_genes, thresh_label = load_degs(deg_file)

        if len(up) + len(down) == 0:
            print(f"  {ct_name}: 0 DEGs ({thresh_label}), skipping")
            continue

        print(f"  {ct_name}: {len(up)} up, {len(down)} down ({thresh_label})")

        ct_outdir = enrichr_dir / ct_name
        ct_outdir.mkdir(exist_ok=True)

        # Run ENRICHR on up and down separately
        if len(up) >= 3:
            res_up = run_enrichr_on_genes(up, f"{ct_name}_Up_Elderly", ct_outdir)
            if len(res_up) > 0:
                res_up["cell_type"] = ct_name
                res_up["direction"] = "Up in Elderly"
                all_enrichr.append(res_up)

        if len(down) >= 3:
            res_down = run_enrichr_on_genes(down, f"{ct_name}_Down_Elderly", ct_outdir)
            if len(res_down) > 0:
                res_down["cell_type"] = ct_name
                res_down["direction"] = "Down in Elderly"
                all_enrichr.append(res_down)

    # -------------------------------------------------------------------------
    # Step 2: Macrophage-specific DEGs (from script 11, relaxed threshold)
    # -------------------------------------------------------------------------
    print("\nStep 2: Loading macrophage-specific DEG results...")
    mac_deg_file = output_dir / "macrophage_degs.csv"
    if mac_deg_file.exists():
        # Try strict FDR first, fall back to nominal p
        up, down, all_genes, _ = load_degs(mac_deg_file, fdr_thresh=FDR_THRESHOLD)
        if len(up) + len(down) == 0:
            print("  0 DEGs at FDR<0.05, using nominal p<0.01")
            up, down, all_genes, thresh_label = load_degs(
                mac_deg_file, use_nominal=True
            )
        else:
            thresh_label = f"FDR<{FDR_THRESHOLD}"

        print(f"  Macrophage-specific: {len(up)} up, {len(down)} down ({thresh_label})")

        mac_outdir = enrichr_dir / "Macrophage_specific"
        mac_outdir.mkdir(exist_ok=True)

        if len(up) >= 3:
            res_up = run_enrichr_on_genes(up, "Macrophage_specific_Up", mac_outdir)
            if len(res_up) > 0:
                res_up["cell_type"] = "Macrophage (specific)"
                res_up["direction"] = "Up in Elderly"
                all_enrichr.append(res_up)

        if len(down) >= 3:
            res_down = run_enrichr_on_genes(
                down, "Macrophage_specific_Down", mac_outdir
            )
            if len(res_down) > 0:
                res_down["cell_type"] = "Macrophage (specific)"
                res_down["direction"] = "Down in Elderly"
                all_enrichr.append(res_down)
    else:
        print("  macrophage_degs.csv not found, skipping")

    # -------------------------------------------------------------------------
    # Step 3: Combine and save results
    # -------------------------------------------------------------------------
    if len(all_enrichr) == 0:
        print("\nNo ENRICHR results. Exiting.")
        sys.exit(0)

    print("\nStep 3: Combining results...")
    enrichr_df = pd.concat(all_enrichr, ignore_index=True)

    # Ensure numeric columns
    for col in ["P-value", "Adjusted P-value", "Odds Ratio", "Combined Score"]:
        if col in enrichr_df.columns:
            enrichr_df[col] = pd.to_numeric(enrichr_df[col], errors="coerce")

    # Save full results
    enrichr_df.to_csv(output_dir / "enrichr_summary.csv", index=False)
    print(f"  Saved enrichr_summary.csv ({len(enrichr_df)} total enrichment results)")

    # Summary: significant pathways per cell type
    sig_enrichr = enrichr_df[enrichr_df["Adjusted P-value"] < 0.05].copy()
    print(f"\n  Significant pathways (adj.p<0.05): {len(sig_enrichr)}")

    if len(sig_enrichr) > 0:
        summary = (
            sig_enrichr.groupby(["cell_type", "direction"])
            .size()
            .reset_index(name="n_pathways")
        )
        print("\n  Pathways per cell type:")
        print(summary.to_string(index=False))

    # -------------------------------------------------------------------------
    # Step 4: Figures
    # -------------------------------------------------------------------------
    print("\nStep 4: Generating figures...")

    import seaborn as sns
    from scipy.cluster.hierarchy import linkage, leaves_list

    def build_enrichr_heatmap(enrichr_df, gene_set_name, title_label, fig_stem,
                              prefix_strip=None, max_pathways=30, min_overlap=2):
        """Build a signed -log10(FDR) heatmap for one ENRICHR library.

        Args:
            enrichr_df: Full ENRICHR results DataFrame.
            gene_set_name: Value of Gene_set column to filter on.
            title_label: Display name for the library (used in title).
            fig_stem: Filename stem (e.g., "fig7c_enrichr_heatmap").
            prefix_strip: Optional prefix to strip from pathway names.
            max_pathways: Max rows to show if too many significant.
            min_overlap: Minimum overlapping genes required (filters spurious
                         single-gene hits common in small pathway databases).
        """
        subset = enrichr_df[
            (enrichr_df["Gene_set"] == gene_set_name)
            & (~enrichr_df["cell_type"].str.contains("specific", case=False))
        ].copy()

        # Parse overlap count from "3/45" format and filter
        if "Overlap" in subset.columns and min_overlap > 1:
            subset["n_overlap"] = subset["Overlap"].str.split("/").str[0].astype(int)
            n_before = len(subset)
            subset = subset[subset["n_overlap"] >= min_overlap]
            n_dropped = n_before - len(subset)
            if n_dropped > 0:
                print(f"    Filtered {n_dropped} terms with <{min_overlap} overlapping genes")

        if len(subset) == 0:
            print(f"  No {title_label} results for heatmap")
            return

        MAX_SCORE = 5  # cap at 5 for balanced color scale (p ~ 1e-5)
        subset["signed_score"] = np.where(
            subset["direction"].str.contains("Up"),
            -np.log10(subset["Adjusted P-value"].clip(lower=1e-30)),
            np.log10(subset["Adjusted P-value"].clip(lower=1e-30)),
        ).clip(-MAX_SCORE, MAX_SCORE)

        # Keep best score per cell type × pathway
        subset["abs_score"] = subset["signed_score"].abs()
        best = subset.sort_values("abs_score", ascending=False).drop_duplicates(
            subset=["cell_type", "Term"], keep="first"
        )

        mat = best.pivot_table(
            index="Term", columns="cell_type", values="signed_score", aggfunc="first"
        ).fillna(0)

        # Clean pathway names
        if prefix_strip:
            mat.index = [
                t.replace(prefix_strip, "").replace("_", " ").title()
                if t.startswith(prefix_strip) else t
                for t in mat.index
            ]
        # Strip trailing identifiers like "Homo sapiens h xxxPathway"
        mat.index = [
            t.split(" Homo sapiens")[0] if " Homo sapiens" in t else t
            for t in mat.index
        ]

        # Select pathways significant in ≥1 cell type
        sig_thresh = -np.log10(0.05)
        sig_mask = (mat.abs() > sig_thresh).any(axis=1)
        mat_sig = mat.loc[sig_mask]

        if len(mat_sig) == 0:
            max_abs = mat.abs().max(axis=1).sort_values(ascending=False)
            mat_sig = mat.loc[max_abs.head(max_pathways).index]

        if len(mat_sig) > max_pathways:
            # Balanced selection: top half by most positive, top half by most negative
            # This prevents one direction from dominating the heatmap
            half = max_pathways // 2
            max_pos = mat_sig.max(axis=1).sort_values(ascending=False)
            max_neg = mat_sig.min(axis=1).sort_values(ascending=True)
            top_pos = set(max_pos[max_pos > 0].head(half).index)
            top_neg = set(max_neg[max_neg < 0].head(half).index)
            keep = list(top_pos | top_neg)
            if len(keep) < max_pathways:
                # Fill remaining with highest absolute score
                remaining = mat_sig.index.difference(keep)
                max_abs = mat_sig.loc[remaining].abs().max(axis=1).sort_values(ascending=False)
                keep.extend(max_abs.head(max_pathways - len(keep)).index.tolist())
            mat_sig = mat_sig.loc[mat_sig.index.isin(keep[:max_pathways])]

        # Sort columns by number of significant pathways
        ct_order = (mat_sig.abs() > sig_thresh).sum(axis=0).sort_values(ascending=False).index
        mat_sig = mat_sig[ct_order]

        # Cluster rows
        if len(mat_sig) > 2:
            row_linkage = linkage(mat_sig.values, method="ward", metric="euclidean")
            mat_sig = mat_sig.iloc[leaves_list(row_linkage)]

        # Plot
        n_pathways = len(mat_sig)
        n_celltypes = len(mat_sig.columns)
        fig_h = max(8, n_pathways * 0.35 + 2)
        fig_w = max(8, n_celltypes * 0.9 + 4)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        vmax = min(MAX_SCORE, mat_sig.abs().values.max())
        sns.heatmap(
            mat_sig, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
            linewidths=0.5, linecolor="white", ax=ax,
            cbar_kws={"label": "Signed −log₁₀(FDR)\n(red = up in Elderly, blue = down)", "shrink": 0.6},
            xticklabels=True, yticklabels=True,
        )

        # Significance markers
        for i, pathway in enumerate(mat_sig.index):
            for j, ct in enumerate(mat_sig.columns):
                val = mat_sig.loc[pathway, ct]
                if abs(val) > -np.log10(0.01):
                    ax.text(j + 0.5, i + 0.5, "**", ha="center", va="center",
                            fontsize=7, fontweight="bold", color="black")
                elif abs(val) > sig_thresh:
                    ax.text(j + 0.5, i + 0.5, "*", ha="center", va="center",
                            fontsize=8, fontweight="bold", color="black")

        ax.set_title(f"ENRICHR {title_label} Pathways: Elderly vs Young\n"
                      "(Pseudobulk DESeq2, * FDR<0.05, ** FDR<0.01)",
                      fontsize=13, pad=12)
        ax.set_xlabel("")
        ax.set_ylabel("")
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=10)
        ax.tick_params(axis="y", rotation=0, labelsize=9)

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.15)
        for ext in ["png", "pdf"]:
            fig.savefig(figures_dir / f"{fig_stem}.{ext}", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved {fig_stem}.png/pdf ({n_pathways} pathways x {n_celltypes} cell types)")

    # 4a: Hallmark heatmap
    print("  Building Hallmark heatmap...")
    build_enrichr_heatmap(
        enrichr_df, "MSigDB_Hallmark_2020", "Hallmark",
        "fig7c_enrichr_heatmap", prefix_strip="HALLMARK_", min_overlap=2,
    )


    # 4b: Macrophage-specific pathways (from macrophage-only pseudobulk, relaxed threshold)
    mac_results = enrichr_df[
        enrichr_df["cell_type"].str.contains("specific", case=False)
    ].copy()
    if len(mac_results) > 0:
        mac_top = mac_results.nsmallest(20, "P-value")
        mac_top["-log10(p)"] = -np.log10(mac_top["P-value"].clip(lower=1e-30))
        mac_top["Term_clean"] = mac_top["Term"].str.replace(
            r"\s*\(.*\)$", "", regex=True
        )

        fig, ax = plt.subplots(figsize=(10, 8))
        colors = [
            "#E41A1C" if "Up" in d else "#377EB8" for d in mac_top["direction"]
        ]
        ax.barh(range(len(mac_top)), mac_top["-log10(p)"], color=colors, alpha=0.8)
        ax.set_yticks(range(len(mac_top)))
        ax.set_yticklabels(
            [row["Term_clean"] for _, row in mac_top.iterrows()], fontsize=9
        )
        ax.set_xlabel("-log10(p-value)", fontsize=12)
        ax.set_title("ENRICHR: Macrophage-Specific Pathways (Elderly vs Young)\n(nominal p<0.01 DEGs)", fontsize=13)
        ax.axvline(-np.log10(0.05), ls="--", color="grey", alpha=0.5, label="p=0.05")
        ax.invert_yaxis()

        from matplotlib.patches import Patch
        ax.legend(handles=[
            Patch(facecolor="#E41A1C", alpha=0.8, label="Up in Elderly"),
            Patch(facecolor="#377EB8", alpha=0.8, label="Down in Elderly"),
        ], loc="lower right")

        plt.tight_layout()
        for ext in ["png", "pdf"]:
            fig.savefig(figures_dir / f"fig7c_enrichr_macrophage.{ext}", dpi=300, bbox_inches="tight")
        plt.close()
        print("  Saved fig7c_enrichr_macrophage.png/pdf")

    # -------------------------------------------------------------------------
    # Step 5: Print key pathway highlights
    # -------------------------------------------------------------------------
    print("\nStep 5: Key pathway highlights...")

    # TNF/TGF pathways (key for manuscript Figure 7C)
    tnf_tgf = sig_enrichr[
        sig_enrichr["Term"].str.contains("TNF|TGF|TNFA|TGFB", case=False, na=False)
    ] if len(sig_enrichr) > 0 else pd.DataFrame()

    if len(tnf_tgf) > 0:
        print("\n  TNF/TGF pathways (significant):")
        for _, row in tnf_tgf.iterrows():
            print(
                f"    {row['cell_type']} | {row['direction']} | "
                f"{row['Term']} | FDR={row['Adjusted P-value']:.2e}"
            )
    else:
        print("  No significant TNF/TGF pathways")

    # Immune/inflammation pathways
    immune_kw = "Inflamm|Immune|Cytokine|Interleukin|Interferon|NF.kB|JAK.STAT"
    immune = sig_enrichr[
        sig_enrichr["Term"].str.contains(immune_kw, case=False, na=False)
    ] if len(sig_enrichr) > 0 else pd.DataFrame()

    if len(immune) > 0:
        print(f"\n  Immune/Inflammation pathways ({len(immune)} significant):")
        for _, row in immune.head(10).iterrows():
            print(
                f"    {row['cell_type']} | {row['direction']} | "
                f"{row['Term']} | FDR={row['Adjusted P-value']:.2e}"
            )

    print("\n=== ENRICHR pathway analysis complete ===")


if __name__ == "__main__":
    main()
