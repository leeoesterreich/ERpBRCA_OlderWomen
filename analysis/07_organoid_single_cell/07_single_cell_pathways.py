#!/usr/bin/env python
# analysis/07_organoid_single_cell/07_single_cell_pathways.py
# Single-cell pathway scoring and distributional comparison across treatments
#
# Uses the single-cell nature of data to compare pathway activity distributions
# between treatment groups. Avoids pseudobulk n=1 limitation.
#
# Inputs:
#   - data/processed/preprocessed.h5ad (from 03_preprocess.py)
#
# Outputs:
#   - data/processed/pathway_scored.h5ad
#   - outputs/single_cell_pathways/pathway_comparisons.csv
#   - figures/single_cell_pathways/violin_*.png/.pdf
#   - figures/single_cell_pathways/ridge_*.png/.pdf
#   - figures/single_cell_pathways/heatmap_effect_sizes.png/.pdf
#   - figures/single_cell_pathways/key_genes_violin.png/.pdf
#   - figures/single_cell_pathways/hypothesis_summary.png/.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings('ignore')

from _config import (
    OUTPUT_DIR, FIGURES_DIR, h5ad_path, check_file_exists, standardize_treatments
)

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_H5AD = h5ad_path("preprocessed.h5ad")
OUTPUT_H5AD = h5ad_path("pathway_scored.h5ad")
RESULTS_DIR = OUTPUT_DIR / "single_cell_pathways"
FIG_DIR = FIGURES_DIR / "single_cell_pathways"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Key comparisons
COMPARISONS = [
    ("E1", "E2", "E1 vs E2"),
    ("E1", "E1+HSD17B7i", "E1 vs E1+HSD17B7i"),
    ("E1+HSD17B7i", "E2+HSD17B7i", "E1+HSD17B7i vs E2+HSD17B7i"),
]

# Gene sets for pathway scoring (from MSigDB Hallmarks and literature)
GENE_SETS = {
    "ESTROGEN_RESPONSE_EARLY": [
        "GREB1", "TFF1", "PGR", "CA12", "XBP1", "CCND1", "MYC", "IGFBP4",
        "PDZK1", "STC2", "NRIP1", "RARA", "KRT19", "ERBB4", "AREG", "EGR3",
        "SIAH2", "PKIB", "HSPB8", "SLC9A3R1", "MYB", "CELSR2", "NPY1R"
    ],
    "ESTROGEN_RESPONSE_LATE": [
        "GREB1", "TFF1", "PGR", "STC2", "PDZK1", "IGFBP4", "CA12", "MYC",
        "CCND1", "XBP1", "NRIP1", "KRT13", "SLC7A2", "MAPT", "ELOVL2",
        "MYBL1", "PRSS23", "FKBP4", "TSKU", "ANXA9"
    ],
    "E2F_TARGETS": [
        "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "RRM1",
        "RRM2", "CDC6", "ORC1", "CDT1", "GINS1", "GINS2", "GINS3", "GINS4",
        "RFC1", "RFC2", "RFC3", "RFC4", "POLE", "POLE2", "POLA1", "PRIM1"
    ],
    "G2M_CHECKPOINT": [
        "CDK1", "CCNB1", "CCNB2", "CDC25C", "PLK1", "AURKA", "AURKB",
        "BUB1", "BUB1B", "MAD2L1", "CENPA", "CENPF", "KIF11", "KIF23",
        "TOP2A", "BIRC5", "CDCA8", "NDC80", "NUF2", "SPC25"
    ],
    "MYC_TARGETS": [
        "MYC", "LDHA", "ENO1", "PKM", "GAPDH", "NPM1", "NCL", "HSPE1",
        "HSPD1", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8",
        "TCP1", "EIF4A1", "EIF4E", "RPS2", "RPL3", "RPL4"
    ],
    "PROLIFERATION": [
        "MKI67", "TOP2A", "PCNA", "MCM2", "CDK1", "CCNB1", "AURKA",
        "PLK1", "BIRC5", "CDC20", "CDCA8", "TPX2", "CENPF", "NUSAP1"
    ],
    "STEROID_BIOSYNTHESIS": [
        "HSD17B7", "HSD17B1", "HSD17B2", "HSD17B4", "CYP19A1", "STS",
        "SULT1E1", "HSD3B1", "HSD3B2", "CYP17A1", "CYP11A1", "STAR"
    ],
    "ER_TARGETS_DIRECT": [
        "GREB1", "TFF1", "PGR", "CCND1", "MYC", "CTSD", "PDZK1", "XBP1",
        "RARA", "CA12", "IGFBP4", "STC2", "AGR2", "FOXA1", "ESR1"
    ]
}

def log_msg(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

def score_pathways(adata):
    """Score cells for each pathway gene set"""
    log_msg("Scoring pathways...")

    for name, genes in GENE_SETS.items():
        # Filter to genes present in data
        genes_present = [g for g in genes if g in adata.var_names]
        if len(genes_present) < 3:
            log_msg(f"  {name}: Only {len(genes_present)} genes found, skipping")
            continue

        sc.tl.score_genes(adata, genes_present, score_name=f"{name}_score")
        log_msg(f"  {name}: {len(genes_present)}/{len(genes)} genes, scored")

    return adata

def compare_distributions(adata, score_col, group1, group2):
    """Compare score distributions between two groups using multiple tests"""
    mask1 = adata.obs["treatment"] == group1
    mask2 = adata.obs["treatment"] == group2

    vals1 = adata.obs.loc[mask1, score_col].dropna()
    vals2 = adata.obs.loc[mask2, score_col].dropna()

    if len(vals1) < 10 or len(vals2) < 10:
        return None

    # Mann-Whitney U test (Wilcoxon rank-sum)
    mw_stat, mw_p = stats.mannwhitneyu(vals1, vals2, alternative='two-sided')

    # Kolmogorov-Smirnov test (tests full distribution shape)
    ks_stat, ks_p = stats.ks_2samp(vals1, vals2)

    # Effect size (Cohen's d)
    pooled_std = np.sqrt((vals1.std()**2 + vals2.std()**2) / 2)
    cohens_d = (vals1.mean() - vals2.mean()) / pooled_std if pooled_std > 0 else 0

    return {
        "group1": group1,
        "group2": group2,
        "n_group1": len(vals1),
        "n_group2": len(vals2),
        "mean_group1": vals1.mean(),
        "mean_group2": vals2.mean(),
        "median_group1": vals1.median(),
        "median_group2": vals2.median(),
        "mannwhitney_stat": mw_stat,
        "mannwhitney_pval": mw_p,
        "ks_stat": ks_stat,
        "ks_pval": ks_p,
        "cohens_d": cohens_d
    }

def run_all_comparisons(adata):
    """Run comparisons for all pathways and treatment pairs"""
    log_msg("Running statistical comparisons...")

    results = []
    score_cols = [c for c in adata.obs.columns if c.endswith("_score") and c != "S_score" and c != "G2M_score"]

    for score_col in score_cols:
        pathway = score_col.replace("_score", "")
        for group1, group2, comp_name in COMPARISONS:
            result = compare_distributions(adata, score_col, group1, group2)
            if result:
                result["pathway"] = pathway
                result["comparison"] = comp_name
                results.append(result)

    df = pd.DataFrame(results)

    # Multiple testing correction (Benjamini-Hochberg)
    if len(df) > 0:
        from scipy.stats import false_discovery_control
        df["mannwhitney_padj"] = false_discovery_control(df["mannwhitney_pval"])
        df["ks_padj"] = false_discovery_control(df["ks_pval"])

    return df

def plot_pathway_violins(adata, pathway_name, score_col):
    """Create violin plot comparing pathway scores across treatments"""
    treatment_order = ["Vehicle", "E1", "E1+fulv", "E1+HSD17B7i", "E2", "E2+fulv", "E2+HSD17B7i"]
    treatment_order = [t for t in treatment_order if t in adata.obs["treatment"].unique()]

    # Color palette - E1 conditions in blues, E2 in reds
    colors = {
        "Vehicle": "#808080",
        "E1": "#3498db",
        "E1+fulv": "#85c1e9",
        "E1+HSD17B7i": "#1a5276",
        "E2": "#e74c3c",
        "E2+fulv": "#f1948a",
        "E2+HSD17B7i": "#922b21"
    }
    palette = [colors.get(t, "#333333") for t in treatment_order]

    fig, ax = plt.subplots(figsize=(10, 5))

    plot_data = adata.obs[["treatment", score_col]].copy()
    plot_data = plot_data[plot_data["treatment"].isin(treatment_order)]

    sns.violinplot(
        data=plot_data, x="treatment", y=score_col,
        order=treatment_order, palette=palette, ax=ax,
        inner="box", cut=0
    )

    ax.set_xlabel("")
    ax.set_ylabel("Pathway Score")
    ax.set_title(f"{pathway_name.replace('_', ' ')} Activity by Treatment")
    plt.xticks(rotation=45, ha="right")

    plt.tight_layout()
    plt.savefig(FIG_DIR / f"violin_{pathway_name}.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIG_DIR / f"violin_{pathway_name}.pdf", bbox_inches="tight")
    plt.close()

def plot_pathway_ridge(adata, pathway_name, score_col):
    """Create ridge plot showing full distributions"""
    treatment_order = ["Vehicle", "E1", "E1+fulv", "E1+HSD17B7i", "E2", "E2+fulv", "E2+HSD17B7i"]
    treatment_order = [t for t in treatment_order if t in adata.obs["treatment"].unique()]

    fig, axes = plt.subplots(len(treatment_order), 1, figsize=(8, len(treatment_order)*1.2),
                             sharex=True)

    colors = {
        "Vehicle": "#808080",
        "E1": "#3498db",
        "E1+fulv": "#85c1e9",
        "E1+HSD17B7i": "#1a5276",
        "E2": "#e74c3c",
        "E2+fulv": "#f1948a",
        "E2+HSD17B7i": "#922b21"
    }

    for i, treatment in enumerate(treatment_order):
        ax = axes[i] if len(treatment_order) > 1 else axes
        mask = adata.obs["treatment"] == treatment
        vals = adata.obs.loc[mask, score_col].dropna()

        ax.fill_between(
            np.linspace(vals.min(), vals.max(), 100),
            0,
            stats.gaussian_kde(vals)(np.linspace(vals.min(), vals.max(), 100)),
            alpha=0.7, color=colors.get(treatment, "#333333")
        )
        ax.axvline(vals.median(), color="black", linestyle="--", alpha=0.5)
        ax.set_ylabel(treatment, rotation=0, ha="right", va="center")
        ax.set_yticks([])

        if i < len(treatment_order) - 1:
            ax.spines["bottom"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    axes[-1].set_xlabel("Pathway Score")
    fig.suptitle(f"{pathway_name.replace('_', ' ')} Distribution", y=1.02)

    plt.tight_layout()
    plt.savefig(FIG_DIR / f"ridge_{pathway_name}.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIG_DIR / f"ridge_{pathway_name}.pdf", bbox_inches="tight")
    plt.close()

def plot_heatmap_summary(results_df):
    """Create heatmap of effect sizes across pathways and comparisons"""
    pivot = results_df.pivot(index="pathway", columns="comparison", values="cohens_d")

    # Reorder comparisons
    comp_order = ["E1 vs E2", "E1 vs E1+HSD17B7i", "E1+HSD17B7i vs E2+HSD17B7i"]
    comp_order = [c for c in comp_order if c in pivot.columns]
    pivot = pivot[comp_order]

    fig, ax = plt.subplots(figsize=(8, max(6, len(pivot)*0.4)))

    sns.heatmap(
        pivot, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
        cbar_kws={"label": "Cohen's d (effect size)"},
        ax=ax, vmin=-1.5, vmax=1.5
    )

    ax.set_xlabel("Comparison")
    ax.set_ylabel("Pathway")
    ax.set_title("Pathway Activity Differences\n(positive = higher in first group)")

    plt.tight_layout()
    plt.savefig(FIG_DIR / "heatmap_effect_sizes.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIG_DIR / "heatmap_effect_sizes.pdf", bbox_inches="tight")
    plt.close()

def plot_key_comparisons(adata):
    """Create focused comparison plots for key hypothesis"""
    log_msg("Creating key comparison plots...")

    # Key genes
    key_genes = ["GREB1", "TFF1", "PGR", "MKI67", "HSD17B7", "ESR1", "CCND1", "MYC"]
    key_genes = [g for g in key_genes if g in adata.var_names]

    treatment_order = ["Vehicle", "E1", "E1+HSD17B7i", "E2", "E2+HSD17B7i"]
    treatment_order = [t for t in treatment_order if t in adata.obs["treatment"].unique()]

    # Dotplot of key genes
    sc.pl.dotplot(
        adata[adata.obs["treatment"].isin(treatment_order)],
        var_names=key_genes,
        groupby="treatment",
        categories_order=treatment_order,
        save="_key_genes.png"
    )

    # Violin plots of key genes
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    colors = {
        "Vehicle": "#808080", "E1": "#3498db", "E1+HSD17B7i": "#1a5276",
        "E2": "#e74c3c", "E2+HSD17B7i": "#922b21"
    }
    palette = [colors[t] for t in treatment_order]

    for i, gene in enumerate(key_genes[:8]):
        if gene in adata.var_names:
            plot_data = adata.obs[["treatment"]].copy()
            plot_data[gene] = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, 'toarray') else adata[:, gene].X.flatten()
            plot_data = plot_data[plot_data["treatment"].isin(treatment_order)]

            sns.violinplot(
                data=plot_data, x="treatment", y=gene,
                order=treatment_order, palette=palette, ax=axes[i],
                inner="box"
            )
            axes[i].set_title(gene)
            axes[i].set_xlabel("")
            axes[i].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "key_genes_violin.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIG_DIR / "key_genes_violin.pdf", bbox_inches="tight")
    plt.close()

def plot_hypothesis_test(adata, results_df):
    """Create summary figure testing the HSD17B7 hypothesis"""
    log_msg("Creating hypothesis test figure...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    treatment_order = ["E1", "E1+HSD17B7i", "E2", "E2+HSD17B7i"]
    treatment_order = [t for t in treatment_order if t in adata.obs["treatment"].unique()]
    colors = {"E1": "#3498db", "E1+HSD17B7i": "#1a5276", "E2": "#e74c3c", "E2+HSD17B7i": "#922b21"}
    palette = [colors[t] for t in treatment_order]

    # Panel A: Estrogen Response Early
    if "ESTROGEN_RESPONSE_EARLY_score" in adata.obs.columns:
        plot_data = adata.obs[["treatment", "ESTROGEN_RESPONSE_EARLY_score"]].copy()
        plot_data = plot_data[plot_data["treatment"].isin(treatment_order)]
        sns.violinplot(data=plot_data, x="treatment", y="ESTROGEN_RESPONSE_EARLY_score",
                      order=treatment_order, palette=palette, ax=axes[0,0], inner="box")
        axes[0,0].set_title("A. Estrogen Response (Early)")
        axes[0,0].set_xlabel("")
        axes[0,0].tick_params(axis='x', rotation=45)

    # Panel B: E2F Targets (cell cycle entry)
    if "E2F_TARGETS_score" in adata.obs.columns:
        plot_data = adata.obs[["treatment", "E2F_TARGETS_score"]].copy()
        plot_data = plot_data[plot_data["treatment"].isin(treatment_order)]
        sns.violinplot(data=plot_data, x="treatment", y="E2F_TARGETS_score",
                      order=treatment_order, palette=palette, ax=axes[0,1], inner="box")
        axes[0,1].set_title("B. E2F Targets (S-phase entry)")
        axes[0,1].set_xlabel("")
        axes[0,1].tick_params(axis='x', rotation=45)

    # Panel C: Proliferation
    if "PROLIFERATION_score" in adata.obs.columns:
        plot_data = adata.obs[["treatment", "PROLIFERATION_score"]].copy()
        plot_data = plot_data[plot_data["treatment"].isin(treatment_order)]
        sns.violinplot(data=plot_data, x="treatment", y="PROLIFERATION_score",
                      order=treatment_order, palette=palette, ax=axes[1,0], inner="box")
        axes[1,0].set_title("C. Proliferation Signature")
        axes[1,0].set_xlabel("")
        axes[1,0].tick_params(axis='x', rotation=45)

    # Panel D: Cell cycle phase proportions
    phase_data = adata.obs[adata.obs["treatment"].isin(treatment_order)]
    phase_props = phase_data.groupby(["treatment", "phase"]).size().unstack(fill_value=0)
    phase_props = phase_props.div(phase_props.sum(axis=1), axis=0)
    phase_props = phase_props.loc[treatment_order]

    phase_props.plot(kind="bar", stacked=True, ax=axes[1,1],
                    color=["#2ecc71", "#f1c40f", "#e74c3c"])
    axes[1,1].set_title("D. Cell Cycle Phase Distribution")
    axes[1,1].set_xlabel("")
    axes[1,1].legend(title="Phase", bbox_to_anchor=(1.02, 1))
    axes[1,1].tick_params(axis='x', rotation=45)

    plt.suptitle("HSD17B7 Inhibitor Effect on Estrogen Signaling\n" +
                 "Hypothesis: Inhibitor reduces E1->E2 conversion, blocking E1 but not E2 effects",
                 fontsize=12, y=1.02)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "hypothesis_summary.png", dpi=150, bbox_inches="tight")
    plt.savefig(FIG_DIR / "hypothesis_summary.pdf", bbox_inches="tight")
    plt.close()

def main():
    log_msg("=" * 60)
    log_msg("Single-cell pathway analysis")
    log_msg("=" * 60)

    # Load data
    log_msg(f"Loading data from {INPUT_H5AD}")
    check_file_exists(INPUT_H5AD, "Preprocessed h5ad")
    adata = sc.read_h5ad(INPUT_H5AD)
    standardize_treatments(adata)
    log_msg(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    log_msg(f"Treatments: {adata.obs['treatment'].value_counts().to_dict()}")

    # Score pathways
    adata = score_pathways(adata)

    # Run statistical comparisons
    results_df = run_all_comparisons(adata)
    results_df.to_csv(RESULTS_DIR / "pathway_comparisons.csv", index=False)
    log_msg(f"Saved comparison results to {RESULTS_DIR / 'pathway_comparisons.csv'}")

    # Print key results
    log_msg("\n" + "=" * 60)
    log_msg("KEY RESULTS (adjusted p < 0.05)")
    log_msg("=" * 60)

    sig_results = results_df[results_df["mannwhitney_padj"] < 0.05].copy()
    sig_results = sig_results.sort_values("cohens_d", key=abs, ascending=False)

    for _, row in sig_results.head(20).iterrows():
        direction = "up" if row["cohens_d"] > 0 else "down"
        log_msg(f"  {row['pathway']:25s} | {row['comparison']:30s} | d={row['cohens_d']:+.2f} {direction} | p={row['mannwhitney_padj']:.2e}")

    # Create visualizations
    log_msg("\nCreating visualizations...")

    # Violin and ridge plots for each pathway
    score_cols = [c for c in adata.obs.columns if c.endswith("_score") and "S_score" not in c and "G2M_score" not in c]
    for score_col in score_cols:
        pathway = score_col.replace("_score", "")
        plot_pathway_violins(adata, pathway, score_col)
        try:
            plot_pathway_ridge(adata, pathway, score_col)
        except:
            pass  # Ridge plots can fail with sparse data

    # Summary heatmap
    if len(results_df) > 0:
        plot_heatmap_summary(results_df)

    # Key gene comparisons
    plot_key_comparisons(adata)

    # Hypothesis test figure
    plot_hypothesis_test(adata, results_df)

    # Save scored data
    adata.write_h5ad(OUTPUT_H5AD)
    log_msg(f"Saved pathway-scored data to {OUTPUT_H5AD}")

    log_msg("\n" + "=" * 60)
    log_msg("Analysis complete!")
    log_msg("=" * 60)

    return adata, results_df

if __name__ == "__main__":
    try:
        adata, results = main()
    except Exception as e:
        print(f"[ERROR] {e}")
        import traceback
        traceback.print_exc()
        raise
