#!/usr/bin/env python
# analysis/07_organoid_single_cell/08_inhibitor_mechanism_exploration.py
# Explore HSD17B7 inhibitor mechanism via signature projection and enzyme profiling
#
# Two competing hypotheses for why the inhibitor affects E2 condition:
#   H1 (Cycling): E2<->E1 interconversion means blocking HSD17B7 traps estrogen
#       as E1 even when exogenous E2 is added
#   H2 (Off-target): Inhibitor hits other HSD17B family members, disrupting
#       steroid metabolism more broadly
#
# Inputs:
#   - data/processed/pathway_scored.h5ad (from 07_single_cell_pathways.py)
#
# Outputs:
#   - data/processed/mechanism_explored.h5ad
#   - outputs/inhibitor_mechanism/estrogen_signature_classification.csv
#   - outputs/inhibitor_mechanism/signature_scores_by_treatment.csv
#   - outputs/inhibitor_mechanism/steroidogenic_gene_expression.csv
#   - outputs/inhibitor_mechanism/dose_response_correlations.csv
#   - outputs/inhibitor_mechanism/subprogram_effects.csv
#   - figures/inhibitor_mechanism/e1_vs_e2_lfc_scatter.png/.pdf
#   - figures/inhibitor_mechanism/signature_projection.png/.pdf
#   - figures/inhibitor_mechanism/heatmap_*.png/.pdf
#   - figures/inhibitor_mechanism/fulv_calibration_kde.png/.pdf
#   - figures/inhibitor_mechanism/dose_response_hsd17b7.png/.pdf
#   - figures/inhibitor_mechanism/subprogram_violins.png/.pdf
#   - figures/inhibitor_mechanism/subprogram_effect_heatmap.png/.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _figure_config import save_fig
import seaborn as sns
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings('ignore')

np.random.seed(42)

from _config import (
    OUTPUT_DIR, FIGURES_DIR, h5ad_path, check_file_exists, standardize_treatments
)

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_H5AD = h5ad_path("pathway_scored.h5ad")
OUTPUT_H5AD = h5ad_path("mechanism_explored.h5ad")
RESULTS_DIR = OUTPUT_DIR / "inhibitor_mechanism"
FIG_DIR = FIGURES_DIR / "inhibitor_mechanism"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Estrogen response sub-programs
SUBPROGRAMS = {
    "ER_PROLIFERATIVE": [
        "MKI67", "TOP2A", "CCND1", "CCNE2", "CDK1", "MCM2", "MCM5",
        "PCNA", "E2F1", "MYC", "MYBL2"
    ],
    "ER_TRANSCRIPTIONAL": [
        "GREB1", "TFF1", "PGR", "PDZK1", "XBP1", "NRIP1", "RARA",
        "FOXA1", "ESR1", "AGR2", "CA12", "STC2"
    ],
    "ER_METABOLIC": [
        "IGFBP4", "PKIB", "SLC7A2", "ELOVL2", "FKBP4", "MAPT",
        "SLC9A3R1", "CELSR2", "PRSS23", "TSKU"
    ],
    "ER_SIGNALING_CROSSTALK": [
        "ERBB4", "AREG", "EGR3", "HSPB8", "NPY1R", "KRT19",
        "SIAH2", "MYB", "CTSD", "ANXA9"
    ]
}

# HSD17B family and steroidogenic enzymes
STEROIDOGENIC_GENES = {
    "HSD17B_FAMILY": [
        "HSD17B1", "HSD17B2", "HSD17B3", "HSD17B4", "HSD17B6", "HSD17B7",
        "HSD17B8", "HSD17B10", "HSD17B11", "HSD17B12", "HSD17B13", "HSD17B14"
    ],
    "OTHER_STEROIDOGENIC": [
        "CYP19A1", "STS", "SULT1E1", "HSD3B1", "HSD3B2",
        "CYP17A1", "CYP11A1", "STAR", "AKR1C1", "AKR1C2", "AKR1C3"
    ],
    "ESTROGEN_RECEPTORS": [
        "ESR1", "ESR2", "GPER1"
    ]
}

def log_msg(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")


# =============================================================================
# MODULE 1: E1-dominant vs E2-dominant signatures
# =============================================================================
def derive_estrogen_signatures(adata):
    """Derive E1-dominant and E2-dominant signatures relative to Vehicle.

    Uses Wilcoxon rank-sum on normalized data to find DE genes.
    Classifies genes as E1-responsive-only, E2-responsive-only, or shared.

    Note: These are equilibrium signatures - both conditions have some
    E1<->E2 interconversion, so they represent net transcriptional states.
    """
    log_msg("MODULE 1: Deriving E1-dominant vs E2-dominant signatures")

    # Work with normalized (log1p) data for DE
    vehicle = adata[adata.obs["treatment"] == "Vehicle"].copy()
    e1 = adata[adata.obs["treatment"] == "E1"].copy()
    e2 = adata[adata.obs["treatment"] == "E2"].copy()

    log_msg(f"  Vehicle: {vehicle.n_obs} cells, E1: {e1.n_obs} cells, E2: {e2.n_obs} cells")

    # Run rank test for E1 vs Vehicle
    e1_v_vehicle = sc.concat([e1, vehicle])
    e1_v_vehicle.obs["group"] = e1_v_vehicle.obs["treatment"].map(
        {"E1": "E1", "Vehicle": "Vehicle"}
    )
    sc.tl.rank_genes_groups(e1_v_vehicle, groupby="group", groups=["E1"],
                            reference="Vehicle", method="wilcoxon")
    e1_de = sc.get.rank_genes_groups_df(e1_v_vehicle, group="E1")
    e1_de = e1_de.rename(columns={"names": "gene", "logfoldchanges": "e1_lfc",
                                   "pvals_adj": "e1_padj", "scores": "e1_score"})

    # Run rank test for E2 vs Vehicle
    e2_v_vehicle = sc.concat([e2, vehicle])
    e2_v_vehicle.obs["group"] = e2_v_vehicle.obs["treatment"].map(
        {"E2": "E2", "Vehicle": "Vehicle"}
    )
    sc.tl.rank_genes_groups(e2_v_vehicle, groupby="group", groups=["E2"],
                            reference="Vehicle", method="wilcoxon")
    e2_de = sc.get.rank_genes_groups_df(e2_v_vehicle, group="E2")
    e2_de = e2_de.rename(columns={"names": "gene", "logfoldchanges": "e2_lfc",
                                   "pvals_adj": "e2_padj", "scores": "e2_score"})

    # Merge results
    de_merged = e1_de[["gene", "e1_lfc", "e1_padj", "e1_score"]].merge(
        e2_de[["gene", "e2_lfc", "e2_padj", "e2_score"]], on="gene", how="outer"
    )

    # Classify genes
    sig_thresh = 0.05
    lfc_thresh = 0.25

    e1_sig = (de_merged["e1_padj"] < sig_thresh) & (de_merged["e1_lfc"].abs() > lfc_thresh)
    e2_sig = (de_merged["e2_padj"] < sig_thresh) & (de_merged["e2_lfc"].abs() > lfc_thresh)

    de_merged["category"] = "not_significant"
    de_merged.loc[e1_sig & ~e2_sig, "category"] = "E1_responsive_only"
    de_merged.loc[~e1_sig & e2_sig, "category"] = "E2_responsive_only"
    de_merged.loc[e1_sig & e2_sig, "category"] = "shared_responsive"

    # Summary
    cat_counts = de_merged["category"].value_counts()
    for cat, count in cat_counts.items():
        log_msg(f"  {cat}: {count} genes")

    # Save
    de_merged.to_csv(RESULTS_DIR / "estrogen_signature_classification.csv", index=False)

    # Extract gene lists for scoring
    e1_up = de_merged[(de_merged["category"] == "E1_responsive_only") &
                       (de_merged["e1_lfc"] > lfc_thresh)]["gene"].tolist()
    e2_up = de_merged[(de_merged["category"] == "E2_responsive_only") &
                       (de_merged["e2_lfc"] > lfc_thresh)]["gene"].tolist()
    shared_up = de_merged[(de_merged["category"] == "shared_responsive") &
                           (de_merged["e1_lfc"] > lfc_thresh) &
                           (de_merged["e2_lfc"] > lfc_thresh)]["gene"].tolist()

    log_msg(f"  E1-dominant upregulated: {len(e1_up)} genes")
    log_msg(f"  E2-dominant upregulated: {len(e2_up)} genes")
    log_msg(f"  Shared upregulated: {len(shared_up)} genes")

    # Save gene lists
    for name, genes in [("e1_dominant_up", e1_up), ("e2_dominant_up", e2_up),
                        ("shared_up", shared_up)]:
        pd.DataFrame({"gene": genes}).to_csv(
            RESULTS_DIR / f"{name}_genes.csv", index=False)

    # Plot: scatter of E1 vs E2 log fold changes
    fig, ax = plt.subplots(figsize=(8, 8))
    colors = {"not_significant": "#cccccc", "E1_responsive_only": "#3498db",
              "E2_responsive_only": "#e74c3c", "shared_responsive": "#9b59b6"}

    for cat in ["not_significant", "shared_responsive", "E1_responsive_only", "E2_responsive_only"]:
        mask = de_merged["category"] == cat
        ax.scatter(de_merged.loc[mask, "e1_lfc"], de_merged.loc[mask, "e2_lfc"],
                  c=colors[cat], s=2, alpha=0.3, label=f"{cat} ({mask.sum()})")

    # Label key genes
    key_genes = ["GREB1", "TFF1", "PGR", "HSD17B7", "ESR1", "MKI67", "CCND1", "MYC"]
    for gene in key_genes:
        row = de_merged[de_merged["gene"] == gene]
        if len(row) > 0:
            ax.annotate(gene, (row["e1_lfc"].values[0], row["e2_lfc"].values[0]),
                       fontsize=10, fontweight="bold")

    ax.axhline(0, color="black", linewidth=0.5, alpha=0.3)
    ax.axvline(0, color="black", linewidth=0.5, alpha=0.3)
    ax.plot([-3, 3], [-3, 3], "k--", alpha=0.3, label="y=x (identical response)")
    ax.set_xlabel("E1 vs Vehicle (log2FC)")
    ax.set_ylabel("E2 vs Vehicle (log2FC)")
    ax.set_title("E1-dominant vs E2-dominant Transcriptional Response\n"
                 "(Note: both conditions have E1<->E2 interconversion)")
    ax.legend(markerscale=5, fontsize=10)
    save_fig(FIG_DIR / "e1_vs_e2_lfc_scatter.png")
    plt.close()

    return de_merged, {"e1_up": e1_up, "e2_up": e2_up, "shared_up": shared_up}


# =============================================================================
# MODULE 2: Signature projection onto inhibitor conditions
# =============================================================================
def project_signatures(adata, gene_lists):
    """Score all conditions with E1-dominant and E2-dominant signatures.

    Key prediction:
      H1 (cycling): E2+HSD17B7i should gain E1-dominant character
      H2 (off-target): Pattern should differ from simple E1 shift
    """
    log_msg("MODULE 2: Projecting signatures onto all conditions")

    # Score with each signature
    for name, genes in gene_lists.items():
        genes_present = [g for g in genes if g in adata.var_names]
        if len(genes_present) >= 3:
            sc.tl.score_genes(adata, genes_present, score_name=f"{name}_signature")
            log_msg(f"  {name}: {len(genes_present)} genes scored")
        else:
            log_msg(f"  {name}: Only {len(genes_present)} genes, skipping")

    # Plot signature scores across all conditions
    treatment_order = ["Vehicle", "E1", "E1+fulv", "E1+HSD17B7i",
                       "E2", "E2+fulv", "E2+HSD17B7i"]
    treatment_order = [t for t in treatment_order if t in adata.obs["treatment"].unique()]

    colors = {
        "Vehicle": "#808080", "E1": "#3498db", "E1+fulv": "#85c1e9",
        "E1+HSD17B7i": "#1a5276", "E2": "#e74c3c", "E2+fulv": "#f1948a",
        "E2+HSD17B7i": "#922b21"
    }

    sig_cols = [c for c in adata.obs.columns if c.endswith("_signature")]
    if len(sig_cols) == 0:
        log_msg("  No signatures to plot")
        return adata

    n_sigs = len(sig_cols)
    fig, axes = plt.subplots(1, n_sigs, figsize=(6 * n_sigs, 5))
    if n_sigs == 1:
        axes = [axes]

    for i, col in enumerate(sig_cols):
        plot_data = adata.obs[["treatment", col]].copy()
        plot_data = plot_data[plot_data["treatment"].isin(treatment_order)]
        palette = [colors.get(t, "#333") for t in treatment_order]

        sns.violinplot(data=plot_data, x="treatment", y=col,
                      order=treatment_order, palette=palette, ax=axes[i],
                      inner="box", cut=0)
        axes[i].set_title(col.replace("_signature", "").replace("_", " ").title() + " Signature")
        axes[i].set_xlabel("")
        axes[i].tick_params(axis="x", rotation=45)

    save_fig(FIG_DIR / "signature_projection.png")
    plt.close()

    # Quantitative summary: mean scores per condition
    summary = adata.obs.groupby("treatment")[sig_cols].mean()
    summary = summary.loc[[t for t in treatment_order if t in summary.index]]
    summary.to_csv(RESULTS_DIR / "signature_scores_by_treatment.csv")
    log_msg("  Signature score summary:")
    for col in sig_cols:
        log_msg(f"    {col}:")
        for t in treatment_order:
            if t in summary.index:
                log_msg(f"      {t:20s}: {summary.loc[t, col]:+.4f}")

    # Key test: E2+HSD17B7i shift toward E1-like?
    # Compare E1-dominant signature: E2 vs E2+HSD17B7i
    if "e1_up_signature" in adata.obs.columns:
        e2_vals = adata.obs.loc[adata.obs["treatment"] == "E2", "e1_up_signature"]
        e2i_vals = adata.obs.loc[adata.obs["treatment"] == "E2+HSD17B7i", "e1_up_signature"]
        if len(e2_vals) > 0 and len(e2i_vals) > 0:
            stat, pval = stats.mannwhitneyu(e2_vals, e2i_vals, alternative="two-sided")
            pooled = np.sqrt((e2_vals.std()**2 + e2i_vals.std()**2) / 2)
            d = (e2_vals.mean() - e2i_vals.mean()) / pooled if pooled > 0 else 0
            log_msg(f"\n  KEY TEST - E1-dominant signature in E2 vs E2+HSD17B7i:")
            log_msg(f"    E2 mean: {e2_vals.mean():.4f}, E2+HSD17B7i mean: {e2i_vals.mean():.4f}")
            log_msg(f"    Cohen's d: {d:.4f}, p={pval:.2e}")
            if e2i_vals.mean() > e2_vals.mean():
                log_msg(f"    -> E2+HSD17B7i shows HIGHER E1-dominant signature (supports cycling hypothesis)")
            else:
                log_msg(f"    -> E2+HSD17B7i shows LOWER E1-dominant signature")

    if "e2_up_signature" in adata.obs.columns:
        e2_vals = adata.obs.loc[adata.obs["treatment"] == "E2", "e2_up_signature"]
        e2i_vals = adata.obs.loc[adata.obs["treatment"] == "E2+HSD17B7i", "e2_up_signature"]
        if len(e2_vals) > 0 and len(e2i_vals) > 0:
            stat, pval = stats.mannwhitneyu(e2_vals, e2i_vals, alternative="two-sided")
            pooled = np.sqrt((e2_vals.std()**2 + e2i_vals.std()**2) / 2)
            d = (e2_vals.mean() - e2i_vals.mean()) / pooled if pooled > 0 else 0
            log_msg(f"\n  KEY TEST - E2-dominant signature in E2 vs E2+HSD17B7i:")
            log_msg(f"    E2 mean: {e2_vals.mean():.4f}, E2+HSD17B7i mean: {e2i_vals.mean():.4f}")
            log_msg(f"    Cohen's d: {d:.4f}, p={pval:.2e}")
            if e2i_vals.mean() < e2_vals.mean():
                log_msg(f"    -> E2+HSD17B7i shows LOWER E2-dominant signature (supports cycling hypothesis)")
            else:
                log_msg(f"    -> E2+HSD17B7i shows HIGHER E2-dominant signature")

    return adata


# =============================================================================
# MODULE 3: HSD17B family profiling
# =============================================================================
def profile_steroidogenic_genes(adata):
    """Profile expression of HSD17B family and steroidogenic enzymes."""
    log_msg("MODULE 3: HSD17B family and steroidogenic enzyme profiling")

    treatment_order = ["Vehicle", "E1", "E1+fulv", "E1+HSD17B7i",
                       "E2", "E2+fulv", "E2+HSD17B7i"]
    treatment_order = [t for t in treatment_order if t in adata.obs["treatment"].unique()]

    all_results = []

    for category, genes in STEROIDOGENIC_GENES.items():
        genes_present = [g for g in genes if g in adata.var_names]
        genes_missing = [g for g in genes if g not in adata.var_names]

        log_msg(f"  {category}: {len(genes_present)}/{len(genes)} genes present")
        if genes_missing:
            log_msg(f"    Missing: {', '.join(genes_missing)}")

        for gene in genes_present:
            for treatment in treatment_order:
                mask = adata.obs["treatment"] == treatment
                expr = adata[mask, gene].X
                if hasattr(expr, 'toarray'):
                    expr = expr.toarray().flatten()
                else:
                    expr = np.array(expr).flatten()

                all_results.append({
                    "category": category,
                    "gene": gene,
                    "treatment": treatment,
                    "mean_expr": expr.mean(),
                    "pct_expressing": (expr > 0).mean() * 100,
                    "median_expr": np.median(expr),
                    "n_cells": mask.sum()
                })

    results_df = pd.DataFrame(all_results)
    results_df.to_csv(RESULTS_DIR / "steroidogenic_gene_expression.csv", index=False)

    # Heatmap of mean expression
    for category, genes in STEROIDOGENIC_GENES.items():
        genes_present = [g for g in genes if g in adata.var_names]
        if len(genes_present) == 0:
            continue

        pivot = results_df[results_df["category"] == category].pivot(
            index="gene", columns="treatment", values="mean_expr"
        )
        pivot = pivot[[t for t in treatment_order if t in pivot.columns]]

        fig, ax = plt.subplots(figsize=(10, max(3, len(genes_present) * 0.5)))
        sns.heatmap(pivot, cmap="YlOrRd", annot=True, fmt=".2f", ax=ax,
                   cbar_kws={"label": "Mean Expression"})
        ax.set_title(f"{category.replace('_', ' ')} Expression Across Treatments")
        save_fig(FIG_DIR / f"heatmap_{category.lower()}.png")
        plt.close()

    # Percent expressing heatmap for HSD17B family
    hsd_genes = [g for g in STEROIDOGENIC_GENES["HSD17B_FAMILY"] if g in adata.var_names]
    if hsd_genes:
        pivot_pct = results_df[results_df["category"] == "HSD17B_FAMILY"].pivot(
            index="gene", columns="treatment", values="pct_expressing"
        )
        pivot_pct = pivot_pct[[t for t in treatment_order if t in pivot_pct.columns]]

        fig, ax = plt.subplots(figsize=(10, max(3, len(hsd_genes) * 0.5)))
        sns.heatmap(pivot_pct, cmap="YlOrRd", annot=True, fmt=".1f", ax=ax,
                   cbar_kws={"label": "% Cells Expressing"})
        ax.set_title("HSD17B Family - % Cells Expressing per Treatment")
        save_fig(FIG_DIR / "heatmap_hsd17b_pct_expressing.png")
        plt.close()

    # Key comparison: expression changes in inhibitor conditions
    log_msg("\n  HSD17B family expression changes with inhibitor:")
    for gene in hsd_genes:
        gene_data = results_df[results_df["gene"] == gene]
        for base, inh in [("E1", "E1+HSD17B7i"), ("E2", "E2+HSD17B7i")]:
            base_row = gene_data[gene_data["treatment"] == base]
            inh_row = gene_data[gene_data["treatment"] == inh]
            if len(base_row) > 0 and len(inh_row) > 0:
                base_expr = base_row["mean_expr"].values[0]
                inh_expr = inh_row["mean_expr"].values[0]
                change = inh_expr - base_expr
                log_msg(f"    {gene:12s} | {base:15s} -> {inh:15s} | "
                       f"{base_expr:.3f} -> {inh_expr:.3f} (delta={change:+.3f})")

    return results_df


# =============================================================================
# MODULE 4: Fulvestrant calibration
# =============================================================================
def fulv_calibration(adata):
    """Use fulvestrant conditions to calibrate blocked vs shifted signaling."""
    log_msg("MODULE 4: Fulvestrant calibration - blocked vs shifted signaling")

    # Get estrogen response scores
    er_scores = [c for c in adata.obs.columns if "ESTROGEN_RESPONSE" in c and c.endswith("_score")]
    if not er_scores:
        log_msg("  No estrogen response scores found, skipping")
        return

    treatment_order = ["Vehicle", "E1", "E1+fulv", "E1+HSD17B7i",
                       "E2", "E2+fulv", "E2+HSD17B7i"]
    treatment_order = [t for t in treatment_order if t in adata.obs["treatment"].unique()]

    # Compare: does E1+HSD17B7i look more like E1+fulv (blocked) or something in between?
    n_plots = min(len(er_scores), 2)
    fig, axes = plt.subplots(1, n_plots, figsize=(7 * n_plots, 5))
    if n_plots == 1:
        axes = [axes]

    for idx, score_col in enumerate(er_scores[:n_plots]):
        pathway_name = score_col.replace("_score", "")

        # Collect distributions
        distributions = {}
        for t in treatment_order:
            mask = adata.obs["treatment"] == t
            distributions[t] = adata.obs.loc[mask, score_col].values

        # KDE overlay
        ax = axes[idx]
        for t, color, ls in [("E1", "#3498db", "-"), ("E1+fulv", "#85c1e9", "--"),
                               ("E1+HSD17B7i", "#1a5276", ":"),
                               ("E2", "#e74c3c", "-"), ("E2+fulv", "#f1948a", "--"),
                               ("E2+HSD17B7i", "#922b21", ":")]:
            if t in distributions and len(distributions[t]) > 0:
                kde = stats.gaussian_kde(distributions[t])
                x = np.linspace(
                    min(distributions[t]) - 0.5,
                    max(distributions[t]) + 0.5, 200
                )
                ax.plot(x, kde(x), color=color, linestyle=ls, label=t, linewidth=2)

        ax.set_title(f"{pathway_name}\nfulv (blocked) vs HSD17B7i (shifted?)")
        ax.set_xlabel("Pathway Score")
        ax.set_ylabel("Density")
        ax.legend(fontsize=10)

    save_fig(FIG_DIR / "fulv_calibration_kde.png")
    plt.close()

    # Quantify: how far is each inhibitor condition from fulv (blocked) vs estrogen (active)?
    log_msg("\n  Inhibitor position between fulv (blocked) and estrogen (active):")
    for score_col in er_scores:
        pathway = score_col.replace("_score", "")
        log_msg(f"  {pathway}:")

        for estrogen, ici, inhib in [("E1", "E1+fulv", "E1+HSD17B7i"),
                                      ("E2", "E2+fulv", "E2+HSD17B7i")]:
            means = {}
            for t in [estrogen, ici, inhib]:
                mask = adata.obs["treatment"] == t
                means[t] = adata.obs.loc[mask, score_col].mean()

            # Position: 0 = same as fulv (fully blocked), 1 = same as estrogen (fully active)
            denom = means[estrogen] - means[ici]
            if abs(denom) > 0.001:
                position = (means[inhib] - means[ici]) / denom
            else:
                position = float('nan')

            log_msg(f"    {inhib:20s}: position = {position:.2f} "
                   f"(0=fulv/blocked, 1=estrogen/active)")
            log_msg(f"      {ici}: {means[ici]:.4f}, {inhib}: {means[inhib]:.4f}, "
                   f"{estrogen}: {means[estrogen]:.4f}")


# =============================================================================
# MODULE 5: Dose-response - HSD17B7 expression vs estrogen response
# =============================================================================
def dose_response_analysis(adata):
    """Within inhibitor conditions, correlate HSD17B7 expression with response."""
    log_msg("MODULE 5: Dose-response - HSD17B7 expression vs estrogen response")

    if "HSD17B7" not in adata.var_names:
        log_msg("  HSD17B7 not in dataset, skipping")
        return

    er_scores = [c for c in adata.obs.columns if "ESTROGEN_RESPONSE" in c and c.endswith("_score")]
    er_scores += [c for c in adata.obs.columns if "ER_TARGETS_DIRECT" in c and c.endswith("_score")]
    er_scores = list(dict.fromkeys(er_scores))  # deduplicate preserving order

    conditions = ["E1+HSD17B7i", "E2+HSD17B7i", "E1", "E2"]
    conditions = [c for c in conditions if c in adata.obs["treatment"].unique()]

    results = []

    n_scores = len(er_scores)
    n_conds = len(conditions)
    if n_scores == 0 or n_conds == 0:
        log_msg("  No scores or conditions to plot")
        return

    fig, axes = plt.subplots(n_conds, n_scores, figsize=(8 * n_scores, 6 * n_conds))
    if n_scores == 1:
        axes = axes.reshape(-1, 1)
    if n_conds == 1:
        axes = axes.reshape(1, -1)

    for i, condition in enumerate(conditions):
        mask = adata.obs["treatment"] == condition
        hsd_expr = adata[mask, "HSD17B7"].X
        if hasattr(hsd_expr, 'toarray'):
            hsd_expr = hsd_expr.toarray().flatten()
        else:
            hsd_expr = np.array(hsd_expr).flatten()

        for j, score_col in enumerate(er_scores):
            er_vals = adata.obs.loc[mask, score_col].values

            # Correlation
            r, p = stats.spearmanr(hsd_expr, er_vals)
            results.append({
                "condition": condition,
                "pathway": score_col,
                "spearman_r": r,
                "spearman_p": p,
                "n_cells": mask.sum(),
                "pct_hsd17b7_positive": (hsd_expr > 0).mean() * 100
            })

            # Plot
            ax = axes[i, j]
            ax.scatter(hsd_expr, er_vals, s=1, alpha=0.1, c="#333333")

            # Bin and show trend
            if (hsd_expr > 0).sum() > 50:
                bins = np.percentile(hsd_expr[hsd_expr > 0],
                                    [0, 25, 50, 75, 100])
                bins = np.unique(np.concatenate([[0], bins]))
                bin_idx = np.digitize(hsd_expr, bins)
                for b in range(1, len(bins)):
                    b_mask = bin_idx == b
                    if b_mask.sum() > 0:
                        ax.scatter(hsd_expr[b_mask].mean(), er_vals[b_mask].mean(),
                                  c="red", s=50, zorder=5, edgecolors="black")

            ax.set_xlabel("HSD17B7 expression", fontsize=16)
            ax.set_ylabel(score_col.replace("_score", ""), fontsize=16)
            ax.set_title(f"{condition}\nr={r:.3f}, p={p:.1e}", fontsize=16)
            ax.tick_params(labelsize=14)

    save_fig(FIG_DIR / "dose_response_hsd17b7.png")
    plt.close()

    results_df = pd.DataFrame(results)
    results_df.to_csv(RESULTS_DIR / "dose_response_correlations.csv", index=False)

    log_msg("  Dose-response correlations:")
    for _, row in results_df.iterrows():
        sig = "***" if row["spearman_p"] < 0.001 else "**" if row["spearman_p"] < 0.01 else "*" if row["spearman_p"] < 0.05 else "ns"
        log_msg(f"    {row['condition']:20s} | {row['pathway']:30s} | "
               f"r={row['spearman_r']:+.3f} {sig} | "
               f"HSD17B7+ cells: {row['pct_hsd17b7_positive']:.1f}%")

    return results_df


# =============================================================================
# MODULE 6: Estrogen response sub-program decomposition
# =============================================================================
def subprogram_decomposition(adata):
    """Break estrogen response into sub-programs and compare effects."""
    log_msg("MODULE 6: Estrogen response sub-program decomposition")

    # Score each sub-program
    for name, genes in SUBPROGRAMS.items():
        genes_present = [g for g in genes if g in adata.var_names]
        if len(genes_present) >= 3:
            sc.tl.score_genes(adata, genes_present, score_name=f"{name}_score")
            log_msg(f"  {name}: {len(genes_present)}/{len(genes)} genes scored")
        else:
            log_msg(f"  {name}: Only {len(genes_present)} genes, skipping")

    treatment_order = ["Vehicle", "E1", "E1+fulv", "E1+HSD17B7i",
                       "E2", "E2+fulv", "E2+HSD17B7i"]
    treatment_order = [t for t in treatment_order if t in adata.obs["treatment"].unique()]

    colors = {
        "Vehicle": "#808080", "E1": "#3498db", "E1+fulv": "#85c1e9",
        "E1+HSD17B7i": "#1a5276", "E2": "#e74c3c", "E2+fulv": "#f1948a",
        "E2+HSD17B7i": "#922b21"
    }

    subprog_scores = [c for c in adata.obs.columns if any(
        c.startswith(sp) for sp in SUBPROGRAMS) and c.endswith("_score")]

    if len(subprog_scores) == 0:
        log_msg("  No sub-program scores computed")
        return adata

    # Violin plots for each sub-program
    n_sp = len(subprog_scores)
    fig, axes = plt.subplots(1, n_sp, figsize=(6 * n_sp, 5))
    if n_sp == 1:
        axes = [axes]

    for i, col in enumerate(subprog_scores):
        plot_data = adata.obs[["treatment", col]].copy()
        plot_data = plot_data[plot_data["treatment"].isin(treatment_order)]
        palette = [colors.get(t, "#333") for t in treatment_order]

        sns.violinplot(data=plot_data, x="treatment", y=col,
                      order=treatment_order, palette=palette, ax=axes[i],
                      inner="box", cut=0)
        axes[i].set_title(col.replace("_score", "").replace("_", " "))
        axes[i].set_xlabel("")
        axes[i].tick_params(axis="x", rotation=45)

    save_fig(FIG_DIR / "subprogram_violins.png")
    plt.close()

    # Effect size heatmap: inhibitor effect on each sub-program
    effect_data = []
    for col in subprog_scores:
        pathway = col.replace("_score", "")
        for group1, group2, label in [
            ("E1", "E1+HSD17B7i", "E1 -> E1+inh"),
            ("E2", "E2+HSD17B7i", "E2 -> E2+inh"),
            ("E1+fulv", "E1+HSD17B7i", "E1+fulv vs E1+inh"),
            ("E2+fulv", "E2+HSD17B7i", "E2+fulv vs E2+inh"),
        ]:
            mask1 = adata.obs["treatment"] == group1
            mask2 = adata.obs["treatment"] == group2
            vals1 = adata.obs.loc[mask1, col].dropna()
            vals2 = adata.obs.loc[mask2, col].dropna()

            if len(vals1) > 10 and len(vals2) > 10:
                pooled_std = np.sqrt((vals1.std()**2 + vals2.std()**2) / 2)
                d = (vals1.mean() - vals2.mean()) / pooled_std if pooled_std > 0 else 0
                _, p = stats.mannwhitneyu(vals1, vals2, alternative="two-sided")
                effect_data.append({
                    "subprogram": pathway,
                    "comparison": label,
                    "cohens_d": d,
                    "pval": p,
                    "mean_group1": vals1.mean(),
                    "mean_group2": vals2.mean()
                })

    effect_df = pd.DataFrame(effect_data)
    effect_df.to_csv(RESULTS_DIR / "subprogram_effects.csv", index=False)

    # Heatmap
    if len(effect_df) > 0:
        pivot = effect_df.pivot(index="subprogram", columns="comparison", values="cohens_d")
        comp_order = ["E1 -> E1+inh", "E2 -> E2+inh", "E1+fulv vs E1+inh", "E2+fulv vs E2+inh"]
        pivot = pivot[[c for c in comp_order if c in pivot.columns]]

        fig, ax = plt.subplots(figsize=(8, max(4, len(pivot) * 0.6)))
        sns.heatmap(pivot, cmap="RdBu_r", center=0, annot=True, fmt=".2f", ax=ax,
                   vmin=-1.0, vmax=1.0,
                   cbar_kws={"label": "Cohen's d (positive = higher in first group)"})
        ax.set_title("Estrogen Sub-program Response to Inhibitor\n"
                    "(Which programs are differentially affected?)")
        save_fig(FIG_DIR / "subprogram_effect_heatmap.png")
        plt.close()

    return adata


# =============================================================================
# MAIN
# =============================================================================
def main():
    log_msg("=" * 60)
    log_msg("HSD17B7 Inhibitor Mechanism Exploration")
    log_msg("=" * 60)

    log_msg(f"Loading data from {INPUT_H5AD}")
    check_file_exists(INPUT_H5AD, "Pathway-scored h5ad")
    adata = sc.read_h5ad(INPUT_H5AD)
    standardize_treatments(adata)
    log_msg(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    log_msg(f"Treatments: {adata.obs['treatment'].value_counts().to_dict()}")

    # Module 1: Derive signatures
    de_results, gene_lists = derive_estrogen_signatures(adata)

    # Module 2: Project signatures
    adata = project_signatures(adata, gene_lists)

    # Module 3: Steroidogenic gene profiling
    steroid_results = profile_steroidogenic_genes(adata)

    # Module 4: Fulvestrant calibration
    fulv_calibration(adata)

    # Module 5: Dose-response
    dose_results = dose_response_analysis(adata)

    # Module 6: Sub-program decomposition
    adata = subprogram_decomposition(adata)

    # Save updated adata
    adata.write_h5ad(OUTPUT_H5AD)
    log_msg(f"\nSaved mechanism-explored data to {OUTPUT_H5AD}")

    log_msg("\n" + "=" * 60)
    log_msg("Analysis complete!")
    log_msg(f"Results: {RESULTS_DIR}/")
    log_msg(f"Figures: {FIG_DIR}/")
    log_msg("=" * 60)

    return adata


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}")
        import traceback
        traceback.print_exc()
        raise
