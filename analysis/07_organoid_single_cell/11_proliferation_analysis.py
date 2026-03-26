#!/usr/bin/env python
# analysis/07_organoid_single_cell/11_proliferation_analysis.py
# Proliferation, quiescence, and cycling heterogeneity analysis
#
# Replaces the naive 3-phase (G1/S/G2M) cell cycle classification with:
# 1. Quiescence/G0 scoring (complements existing PROLIFERATION_score from script 07)
# 2. Binary cycling classification via GMM on PROLIFERATION_score
# 3. Statistical comparisons with effect sizes (Wilcoxon + Fisher's exact)
# 4. Per-cluster heterogeneity analysis
# 5. Proliferation-estrogen response coupling analysis
#
# Inputs:
#   - data/processed/pathway_scored.h5ad (from 07_single_cell_pathways.py)
#
# Outputs:
#   - data/processed/proliferation_analyzed.h5ad
#   - outputs/proliferation/score_comparisons.csv
#   - outputs/proliferation/cycling_fraction_comparisons.csv
#   - outputs/proliferation/per_cluster_cycling.csv
#   - outputs/proliferation/proliferation_estrogen_correlation.csv
#   - figures/proliferation/panel1_score_violins.png/.pdf
#   - figures/proliferation/panel2_cycling_fraction.png/.pdf
#   - figures/proliferation/panel3_key_markers.png/.pdf
#   - figures/proliferation/panel4_score_ridges.png/.pdf
#   - figures/proliferation/panel5a_umap_scores.png/.pdf
#   - figures/proliferation/panel5b_umap_cycling_by_treatment.png/.pdf
#   - figures/proliferation/panel6_cluster_cycling_heatmap.png/.pdf
#   - figures/proliferation/panel7_prolif_er_coupling.png/.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _figure_config import save_fig
import seaborn as sns
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu, fisher_exact, false_discovery_control, spearmanr
from sklearn.mixture import GaussianMixture
from scipy.sparse import issparse
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings('ignore')

from _config import (
    OUTPUT_DIR, FIGURES_DIR, h5ad_path, check_file_exists, standardize_treatments,
    TREATMENT_COLORS, TREATMENT_ORDER,
)

# =============================================================================
# CONFIGURATION
# =============================================================================
INPUT_H5AD = h5ad_path("pathway_scored.h5ad")
OUTPUT_H5AD = h5ad_path("proliferation_analyzed.h5ad")
RESULTS_DIR = OUTPUT_DIR / "proliferation"
FIG_DIR = FIGURES_DIR / "proliferation"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Quiescence/G0 gene set (RB1/HES1 excluded per spec)
QUIESCENCE_GENES = ['CDKN1A', 'CDKN1B', 'BTG1', 'BTG2', 'TOB1', 'GAS1', 'CDKN2A', 'CCNG2']

# Key marker genes for individual violin plots
KEY_MARKERS = ['MKI67', 'PCNA', 'TOP2A', 'CDKN1A', 'CDKN1B']

# Statistical comparisons
COMPARISONS = [
    ('Vehicle', 'E1', 'Vehicle vs E1'),
    ('Vehicle', 'E2', 'Vehicle vs E2'),
    ('E1', 'E2', 'E1 vs E2'),
    ('E1', 'E1+fulv', 'E1 vs E1+fulv'),
    ('E2', 'E2+fulv', 'E2 vs E2+fulv'),
    ('E1', 'E1+HSD17B7i', 'E1 vs E1+HSD17B7i'),
    ('E1+HSD17B7i', 'E2+HSD17B7i', 'E1+HSD17B7i vs E2+HSD17B7i'),
]


def log_msg(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")


# =============================================================================
# SCORING FUNCTIONS
# =============================================================================

def score_quiescence(adata):
    """Score cells for quiescence/G0 using arrest markers.

    Uses: CDKN1A (p21), CDKN1B (p27), BTG1, BTG2, TOB1, GAS1, CDKN2A (p16), CCNG2.
    RB1 excluded (mutation-prone in breast cancer), HES1 excluded (Notch confound).

    Parameters
    ----------
    adata : AnnData
        Must have gene expression data with at least some quiescence genes.

    Returns
    -------
    adata : AnnData
        Adds 'quiescence_score' to adata.obs.
    """
    genes_present = [g for g in QUIESCENCE_GENES if g in adata.var_names]
    genes_missing = [g for g in QUIESCENCE_GENES if g not in adata.var_names]

    log_msg(f"Scoring quiescence: {len(genes_present)}/{len(QUIESCENCE_GENES)} genes found")
    if genes_missing:
        log_msg(f"  Missing: {genes_missing}")

    if len(genes_present) < 3:
        raise ValueError(f"Too few quiescence genes found ({len(genes_present)}). Need at least 3.")

    sc.tl.score_genes(adata, genes_present, score_name='quiescence_score')

    log_msg(f"  Quiescence score range: [{adata.obs['quiescence_score'].min():.3f}, {adata.obs['quiescence_score'].max():.3f}]")
    log_msg(f"  Quiescence score mean: {adata.obs['quiescence_score'].mean():.3f}")

    return adata


def classify_cycling(adata):
    """Classify cells as cycling vs non-cycling using GMM on PROLIFERATION_score.

    Fits a 2-component Gaussian mixture model globally (all cells, all conditions).
    Validates bimodality via BIC comparison (1 vs 2 components).
    Falls back to 75th percentile if GMM fails or BIC favors 1 component.

    Parameters
    ----------
    adata : AnnData
        Must have 'PROLIFERATION_score' in adata.obs.

    Returns
    -------
    adata : AnnData
        Adds 'cycling_status' ('cycling'/'non-cycling'), 'cycling_phase', and
        'cycling_threshold' to adata.obs/uns.
    """
    scores = adata.obs['PROLIFERATION_score'].values.reshape(-1, 1)

    # Fit 1-component and 2-component GMM
    gmm1 = GaussianMixture(n_components=1, random_state=42).fit(scores)
    gmm2 = GaussianMixture(n_components=2, random_state=42).fit(scores)

    bic1 = gmm1.bic(scores)
    bic2 = gmm2.bic(scores)

    log_msg(f"GMM BIC: 1-component={bic1:.1f}, 2-component={bic2:.1f}")

    if bic2 < bic1 and gmm2.converged_:
        # Use GMM threshold: midpoint between the two component means
        means = gmm2.means_.flatten()
        threshold = np.mean(means)
        method = 'GMM'
        log_msg(f"  GMM converged. Component means: {means[0]:.3f}, {means[1]:.3f}")
    else:
        # Fallback: 75th percentile
        threshold = np.percentile(adata.obs['PROLIFERATION_score'], 75)
        method = '75th percentile'
        log_msg(f"  GMM did not improve fit. Using 75th percentile fallback.")

    log_msg(f"  Cycling threshold: {threshold:.3f} (method: {method})")

    # Classify
    adata.obs['cycling_status'] = np.where(
        adata.obs['PROLIFERATION_score'] >= threshold, 'cycling', 'non-cycling'
    )

    # Store threshold info
    adata.uns['cycling_threshold'] = threshold
    adata.uns['cycling_method'] = method

    # Log summary
    cycling_counts = adata.obs['cycling_status'].value_counts()
    total = len(adata)
    for status in ['cycling', 'non-cycling']:
        n = cycling_counts.get(status, 0)
        log_msg(f"  {status}: {n} cells ({100*n/total:.1f}%)")

    # Sub-classify cycling cells into S vs G2M using existing Tirosh scores
    cycling_mask = adata.obs['cycling_status'] == 'cycling'
    adata.obs['cycling_phase'] = 'non-cycling'
    if 'S_score' in adata.obs.columns and 'G2M_score' in adata.obs.columns:
        adata.obs.loc[cycling_mask, 'cycling_phase'] = np.where(
            adata.obs.loc[cycling_mask, 'S_score'] > adata.obs.loc[cycling_mask, 'G2M_score'],
            'S', 'G2M'
        )
        phase_counts = adata.obs.loc[cycling_mask, 'cycling_phase'].value_counts()
        log_msg(f"  Cycling sub-classification: S={phase_counts.get('S', 0)}, G2M={phase_counts.get('G2M', 0)}")

    # Per-treatment breakdown
    log_msg("  Per-treatment cycling fractions:")
    for treatment in TREATMENT_ORDER:
        mask = adata.obs['treatment'] == treatment
        if mask.sum() == 0:
            continue
        n_cycling = (adata.obs.loc[mask, 'cycling_status'] == 'cycling').sum()
        n_total = mask.sum()
        log_msg(f"    {treatment}: {n_cycling}/{n_total} ({100*n_cycling/n_total:.1f}%)")

    return adata


# =============================================================================
# STATISTICAL COMPARISON FUNCTIONS
# =============================================================================

def compare_scores(adata, score_col, group1, group2):
    """Compare a continuous score between two treatment groups.

    Uses Wilcoxon rank-sum test with rank-biserial correlation as effect size.

    Parameters
    ----------
    adata : AnnData
    score_col : str
        Column name in adata.obs
    group1, group2 : str
        Treatment names

    Returns
    -------
    dict with test results, or None if insufficient data
    """
    vals1 = adata.obs.loc[adata.obs['treatment'] == group1, score_col].dropna()
    vals2 = adata.obs.loc[adata.obs['treatment'] == group2, score_col].dropna()

    if len(vals1) < 10 or len(vals2) < 10:
        return None

    stat, pval = mannwhitneyu(vals1, vals2, alternative='two-sided')

    # Rank-biserial correlation: positive r = group1 > group2
    n1, n2 = len(vals1), len(vals2)
    rank_biserial = (2 * stat) / (n1 * n2) - 1

    return {
        'score': score_col,
        'group1': group1,
        'group2': group2,
        'n1': n1,
        'n2': n2,
        'median1': vals1.median(),
        'median2': vals2.median(),
        'mean1': vals1.mean(),
        'mean2': vals2.mean(),
        'U_stat': stat,
        'pval': pval,
        'rank_biserial': rank_biserial,
    }


def compare_cycling_fraction(adata, group1, group2):
    """Compare cycling fractions between two treatment groups.

    Uses Fisher's exact test with odds ratio as effect size.

    Parameters
    ----------
    adata : AnnData
        Must have 'cycling_status' and 'treatment' in obs.
    group1, group2 : str
        Treatment names

    Returns
    -------
    dict with test results, or None if insufficient data
    """
    for g in [group1, group2]:
        if g not in adata.obs['treatment'].values:
            return None

    mask1 = adata.obs['treatment'] == group1
    mask2 = adata.obs['treatment'] == group2

    cyc1 = (adata.obs.loc[mask1, 'cycling_status'] == 'cycling').sum()
    noncyc1 = (adata.obs.loc[mask1, 'cycling_status'] == 'non-cycling').sum()
    cyc2 = (adata.obs.loc[mask2, 'cycling_status'] == 'cycling').sum()
    noncyc2 = (adata.obs.loc[mask2, 'cycling_status'] == 'non-cycling').sum()

    table = np.array([[cyc1, noncyc1], [cyc2, noncyc2]])
    odds_ratio, pval = fisher_exact(table)

    return {
        'group1': group1,
        'group2': group2,
        'cycling_1': cyc1,
        'noncycling_1': noncyc1,
        'frac_cycling_1': cyc1 / (cyc1 + noncyc1),
        'cycling_2': cyc2,
        'noncycling_2': noncyc2,
        'frac_cycling_2': cyc2 / (cyc2 + noncyc2),
        'odds_ratio': odds_ratio,
        'pval': pval,
    }


def run_all_statistics(adata):
    """Run all statistical comparisons for proliferation, quiescence, and cycling fraction.

    Returns
    -------
    score_results : DataFrame
        Wilcoxon results for proliferation_score and quiescence_score
    cycling_results : DataFrame
        Fisher's exact results for cycling fraction
    """
    log_msg("\nRunning statistical comparisons...")

    # Score comparisons (proliferation + quiescence)
    score_rows = []
    for score_col in ['PROLIFERATION_score', 'quiescence_score']:
        for g1, g2, comp_name in COMPARISONS:
            result = compare_scores(adata, score_col, g1, g2)
            if result:
                result['comparison'] = comp_name
                score_rows.append(result)

    score_df = pd.DataFrame(score_rows)

    # BH correction per metric
    if len(score_df) > 0:
        for score_col in score_df['score'].unique():
            mask = score_df['score'] == score_col
            score_df.loc[mask, 'padj'] = false_discovery_control(score_df.loc[mask, 'pval'])

    # Cycling fraction comparisons
    cycling_rows = []
    for g1, g2, comp_name in COMPARISONS:
        result = compare_cycling_fraction(adata, g1, g2)
        if result:
            result['comparison'] = comp_name
            cycling_rows.append(result)

    cycling_df = pd.DataFrame(cycling_rows)
    if len(cycling_df) > 0:
        cycling_df['padj'] = false_discovery_control(cycling_df['pval'])

    # Log key results
    log_msg("\nScore comparison results (padj < 0.05):")
    if len(score_df) > 0:
        sig = score_df[score_df['padj'] < 0.05]
        for _, row in sig.iterrows():
            log_msg(f"  {row['score']:25s} | {row['comparison']:30s} | r={row['rank_biserial']:+.3f} | padj={row['padj']:.2e}")

    log_msg("\nCycling fraction results (padj < 0.05):")
    if len(cycling_df) > 0:
        sig = cycling_df[cycling_df['padj'] < 0.05]
        for _, row in sig.iterrows():
            log_msg(f"  {row['comparison']:30s} | OR={row['odds_ratio']:.2f} | {row['frac_cycling_1']:.1%} vs {row['frac_cycling_2']:.1%} | padj={row['padj']:.2e}")

    return score_df, cycling_df


# =============================================================================
# VISUALIZATION FUNCTIONS -- CORE (Panels 1-4)
# =============================================================================

def plot_score_violins(adata):
    """Panel 1: Violin/box plots of proliferation and quiescence scores by treatment."""
    log_msg("Creating Panel 1: Score violin plots...")

    available = [t for t in TREATMENT_ORDER if t in adata.obs['treatment'].unique()]
    palette = [TREATMENT_COLORS[t] for t in available]

    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    for ax, score_col, title in zip(
        axes,
        ['PROLIFERATION_score', 'quiescence_score'],
        ['Proliferation Score', 'Quiescence (G0) Score']
    ):
        plot_data = adata.obs[['treatment', score_col]].copy()
        plot_data = plot_data[plot_data['treatment'].isin(available)]

        sns.violinplot(
            data=plot_data, x='treatment', y=score_col,
            order=available, palette=palette, ax=ax,
            inner='box', cut=0
        )
        ax.set_xlabel('')
        ax.set_ylabel(title, fontsize=11)
        ax.set_title(title + ' by Treatment', fontsize=12)
        ax.tick_params(axis='x', rotation=30)

        # Add cell counts
        for i, t in enumerate(available):
            n = (adata.obs['treatment'] == t).sum()
            ax.text(i, ax.get_ylim()[0] - 0.02 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                    f'n={n}', ha='center', va='top', fontsize=10, color='gray')

    save_fig(FIG_DIR / "panel1_score_violins.png")
    plt.close()
    log_msg(f"  Saved: panel1_score_violins.png/.pdf")


def plot_cycling_fraction(adata):
    """Panel 2: Stacked bar plot of cycling vs non-cycling per treatment."""
    log_msg("Creating Panel 2: Cycling fraction bar plot...")

    available = [t for t in TREATMENT_ORDER if t in adata.obs['treatment'].unique()]

    counts = adata.obs.groupby(['treatment', 'cycling_status']).size().unstack(fill_value=0)
    props = counts.div(counts.sum(axis=1), axis=0)
    props = props.loc[available]
    counts = counts.loc[available]

    fig, ax = plt.subplots(figsize=(10, 5))

    x = np.arange(len(available))
    cycling_vals = props.get('cycling', pd.Series(0, index=available)).values
    noncycling_vals = props.get('non-cycling', pd.Series(0, index=available)).values

    ax.bar(x, noncycling_vals, label='Non-cycling', color='#3498db', edgecolor='white', linewidth=0.5)
    ax.bar(x, cycling_vals, bottom=noncycling_vals, label='Cycling', color='#e74c3c', edgecolor='white', linewidth=0.5)

    ax.set_xlabel('')
    ax.set_ylabel('Proportion', fontsize=11)
    ax.set_title('Cycling vs Non-Cycling Cells by Treatment', fontsize=12, pad=14)
    ax.set_xticks(x)
    ax.set_xticklabels(available, rotation=30, ha='right')
    ax.set_ylim(0, 1)
    ax.legend(loc='upper right')

    # Annotate cell counts
    for i, t in enumerate(available):
        n_total = counts.loc[t].sum()
        n_cyc = counts.loc[t].get('cycling', 0)
        ax.text(i, 0.98, f'n={n_total}\n({n_cyc} cyc)', ha='center', va='top', fontsize=10, color='gray')

    save_fig(FIG_DIR / "panel2_cycling_fraction.png")
    plt.close()
    log_msg(f"  Saved: panel2_cycling_fraction.png/.pdf")


def plot_key_markers(adata):
    """Panel 3: Individual gene violin plots for key proliferation/arrest markers."""
    log_msg("Creating Panel 3: Key marker violin plots...")

    available_treatments = [t for t in TREATMENT_ORDER if t in adata.obs['treatment'].unique()]
    palette = [TREATMENT_COLORS[t] for t in available_treatments]

    markers_present = [g for g in KEY_MARKERS if g in adata.var_names]
    markers_missing = [g for g in KEY_MARKERS if g not in adata.var_names]
    if markers_missing:
        log_msg(f"  Missing markers: {markers_missing}")

    n_markers = len(markers_present)
    if n_markers == 0:
        log_msg("  No key markers found in data, skipping Panel 3")
        return

    ncols = min(n_markers, 3)
    nrows = (n_markers + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    if n_markers == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for i, gene in enumerate(markers_present):
        ax = axes[i]
        gene_vals = adata[:, gene].X.toarray().flatten() if issparse(adata[:, gene].X) else adata[:, gene].X.flatten()
        plot_data = pd.DataFrame({
            'treatment': adata.obs['treatment'].values,
            gene: gene_vals
        })
        plot_data = plot_data[plot_data['treatment'].isin(available_treatments)]

        sns.violinplot(
            data=plot_data, x='treatment', y=gene,
            order=available_treatments, palette=palette, ax=ax,
            inner='box', cut=0
        )
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('Expression')
        ax.tick_params(axis='x', rotation=30)

    # Hide unused axes
    for j in range(n_markers, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle('Key Proliferation & Arrest Markers by Treatment', fontsize=13, y=1.02)
    save_fig(FIG_DIR / "panel3_key_markers.png")
    plt.close()
    log_msg(f"  Saved: panel3_key_markers.png/.pdf")


def plot_score_ridges(adata):
    """Panel 4: Ridge/density plots of proliferation score per treatment, overlaid."""
    log_msg("Creating Panel 4: Proliferation score ridge plots...")

    available = [t for t in TREATMENT_ORDER if t in adata.obs['treatment'].unique()]

    fig, axes = plt.subplots(len(available), 1, figsize=(8, len(available) * 1.4), sharex=True)

    for i, treatment in enumerate(available):
        ax = axes[i] if len(available) > 1 else axes
        mask = adata.obs['treatment'] == treatment
        vals = adata.obs.loc[mask, 'PROLIFERATION_score'].dropna()

        # KDE
        x_range = np.linspace(
            adata.obs['PROLIFERATION_score'].min() - 0.1,
            adata.obs['PROLIFERATION_score'].max() + 0.1,
            200
        )
        kde = stats.gaussian_kde(vals)
        density = kde(x_range)

        ax.fill_between(x_range, 0, density, alpha=0.7, color=TREATMENT_COLORS.get(treatment, '#333'))
        ax.axvline(vals.median(), color='black', linestyle='--', alpha=0.5, linewidth=1)

        # Mark cycling threshold
        if 'cycling_threshold' in adata.uns:
            ax.axvline(adata.uns['cycling_threshold'], color='red', linestyle=':', alpha=0.5, linewidth=1)

        ax.set_ylabel(treatment, rotation=0, ha='right', va='center', fontsize=10)
        ax.set_yticks([])

        for spine in ['top', 'right', 'left']:
            ax.spines[spine].set_visible(False)
        if i < len(available) - 1:
            ax.spines['bottom'].set_visible(False)

    axes[-1].set_xlabel('Proliferation Score')
    fig.suptitle('Proliferation Score Distributions by Treatment\n(dashed black = median, dotted red = cycling threshold)',
                 fontsize=11, y=1.02)

    save_fig(FIG_DIR / "panel4_score_ridges.png")
    plt.close()
    log_msg(f"  Saved: panel4_score_ridges.png/.pdf")


# =============================================================================
# VISUALIZATION FUNCTIONS -- HETEROGENEITY (Panels 5-7)
# =============================================================================

def plot_umap_overlays(adata):
    """Panel 5: UMAP overlays of proliferation, quiescence, MKI67, and cycling status."""
    log_msg("Creating Panel 5: UMAP overlays...")

    if 'X_umap' not in adata.obsm:
        log_msg("  No UMAP coordinates found, skipping Panel 5")
        return

    # 5a: Score overlays (proliferation, quiescence, MKI67)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    umap = adata.obsm['X_umap']

    # Proliferation score
    sc_plot = axes[0].scatter(umap[:, 0], umap[:, 1], c=adata.obs['PROLIFERATION_score'],
                               cmap='YlOrRd', s=1, alpha=0.5, rasterized=True)
    plt.colorbar(sc_plot, ax=axes[0], shrink=0.8)
    axes[0].set_title('Proliferation Score')
    axes[0].set_xlabel('UMAP1')
    axes[0].set_ylabel('UMAP2')

    # Quiescence score
    sc_plot = axes[1].scatter(umap[:, 0], umap[:, 1], c=adata.obs['quiescence_score'],
                               cmap='YlGnBu', s=1, alpha=0.5, rasterized=True)
    plt.colorbar(sc_plot, ax=axes[1], shrink=0.8)
    axes[1].set_title('Quiescence (G0) Score')
    axes[1].set_xlabel('UMAP1')
    axes[1].set_ylabel('')

    # MKI67 expression
    if 'MKI67' in adata.var_names:
        mki67 = adata[:, 'MKI67'].X.toarray().flatten() if issparse(adata[:, 'MKI67'].X) else adata[:, 'MKI67'].X.flatten()
        sc_plot = axes[2].scatter(umap[:, 0], umap[:, 1], c=mki67,
                                   cmap='Purples', s=1, alpha=0.5, rasterized=True)
        plt.colorbar(sc_plot, ax=axes[2], shrink=0.8)
        axes[2].set_title('MKI67 Expression')
    else:
        axes[2].text(0.5, 0.5, 'MKI67 not found', transform=axes[2].transAxes, ha='center')
    axes[2].set_xlabel('UMAP1')
    axes[2].set_ylabel('')

    save_fig(FIG_DIR / "panel5a_umap_scores.png")
    plt.close()

    # 5b: Treatment-split UMAPs colored by cycling status
    available = [t for t in TREATMENT_ORDER if t in adata.obs['treatment'].unique()]
    ncols = 4
    nrows = (len(available) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = axes.flatten()

    cycling_colors = {'cycling': '#e74c3c', 'non-cycling': '#3498db'}

    for i, treatment in enumerate(available):
        ax = axes[i]
        mask = adata.obs['treatment'] == treatment

        # Background: all cells in light gray
        ax.scatter(umap[:, 0], umap[:, 1], c='#eeeeee', s=0.5, alpha=0.3, rasterized=True)

        # Foreground: this treatment colored by cycling status
        for status, color in cycling_colors.items():
            sub_mask = mask & (adata.obs['cycling_status'] == status)
            if sub_mask.sum() > 0:
                ax.scatter(umap[sub_mask, 0], umap[sub_mask, 1], c=color, s=1, alpha=0.5,
                          label=status, rasterized=True)

        n_cyc = (mask & (adata.obs['cycling_status'] == 'cycling')).sum()
        n_total = mask.sum()
        ax.set_title(f"{treatment}\n({n_cyc}/{n_total} cycling, {100*n_cyc/n_total:.0f}%)", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])

    # Hide unused axes
    for j in range(len(available), len(axes)):
        axes[j].set_visible(False)

    # Add shared legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#e74c3c', label='Cycling'),
                       Patch(facecolor='#3498db', label='Non-cycling')]
    fig.legend(handles=legend_elements, loc='lower right', fontsize=10)

    plt.suptitle('Cycling Status by Treatment (UMAP)', fontsize=12, y=1.01)
    save_fig(FIG_DIR / "panel5b_umap_cycling_by_treatment.png")
    plt.close()
    log_msg(f"  Saved: panel5a_umap_scores.png/.pdf")
    log_msg(f"  Saved: panel5b_umap_cycling_by_treatment.png/.pdf")


def plot_cluster_cycling_heatmap(adata):
    """Panel 6: Heatmap of cycling fraction per Leiden cluster per treatment.

    Identifies resistant subpopulations that keep cycling under fulvestrant.
    """
    log_msg("Creating Panel 6: Per-cluster cycling heatmap...")

    if 'leiden' not in adata.obs.columns:
        log_msg("  No leiden clusters found, skipping Panel 6")
        return None

    available = [t for t in TREATMENT_ORDER if t in adata.obs['treatment'].unique()]

    # Calculate cycling fraction per cluster per treatment
    cluster_cycling = []
    clusters = sorted(adata.obs['leiden'].unique(), key=lambda x: int(x) if x.isdigit() else x)

    for cluster in clusters:
        for treatment in available:
            mask = (adata.obs['leiden'] == cluster) & (adata.obs['treatment'] == treatment)
            n_total = mask.sum()
            if n_total == 0:
                continue
            n_cycling = (adata.obs.loc[mask, 'cycling_status'] == 'cycling').sum()
            cluster_cycling.append({
                'cluster': cluster,
                'treatment': treatment,
                'n_total': n_total,
                'n_cycling': n_cycling,
                'frac_cycling': n_cycling / n_total,
            })

    cluster_df = pd.DataFrame(cluster_cycling)

    # Save CSV
    cluster_df.to_csv(RESULTS_DIR / "per_cluster_cycling.csv", index=False)

    # Pivot for heatmap — fill missing combos with 0 (no cells = 0 cycling)
    pivot = cluster_df.pivot(index='cluster', columns='treatment', values='frac_cycling')
    pivot = pivot.reindex(columns=available).fillna(0)

    fig, ax = plt.subplots(figsize=(10, max(4, len(pivot) * 0.5)))

    annot_labels = pivot.copy().applymap(lambda x: f'{x:.2f}')

    sns.heatmap(
        pivot, cmap='RdYlBu_r', vmin=0, vmax=1, annot=annot_labels, fmt='s',
        cbar_kws={'label': 'Fraction Cycling'}, ax=ax,
        linewidths=0.5, linecolor='white'
    )
    ax.set_xlabel('Treatment')
    ax.set_ylabel('Leiden Cluster')
    ax.set_title('Cycling Fraction per Cluster per Treatment')

    save_fig(FIG_DIR / "panel6_cluster_cycling_heatmap.png")
    plt.close()
    log_msg(f"  Saved: panel6_cluster_cycling_heatmap.png/.pdf")
    log_msg(f"  Saved: {RESULTS_DIR / 'per_cluster_cycling.csv'}")

    return cluster_df


def plot_proliferation_estrogen_coupling(adata):
    """Panel 7: Scatter of proliferation vs estrogen response score, colored by treatment.

    Reports per-treatment Spearman correlation.
    Tests whether HSD17B7i decouples proliferation from estrogen signaling.
    """
    log_msg("Creating Panel 7: Proliferation-estrogen coupling...")

    # Find best estrogen response score column
    er_score_col = None
    for candidate in ['ESTROGEN_RESPONSE_EARLY_score', 'ER_TARGETS_DIRECT_score', 'ESTROGEN_RESPONSE_LATE_score']:
        if candidate in adata.obs.columns:
            er_score_col = candidate
            break

    if er_score_col is None:
        log_msg("  No estrogen response score found, skipping Panel 7")
        return None

    log_msg(f"  Using {er_score_col} for estrogen response axis")

    available = [t for t in TREATMENT_ORDER if t in adata.obs['treatment'].unique()]

    # Calculate per-treatment Spearman correlations
    corr_results = []
    for treatment in available:
        mask = adata.obs['treatment'] == treatment
        prolif = adata.obs.loc[mask, 'PROLIFERATION_score']
        er = adata.obs.loc[mask, er_score_col]
        rho, pval = spearmanr(prolif, er)
        corr_results.append({
            'treatment': treatment,
            'spearman_rho': rho,
            'pval': pval,
            'n': mask.sum(),
        })
        log_msg(f"    {treatment}: rho={rho:.3f}, p={pval:.2e}, n={mask.sum()}")

    corr_df = pd.DataFrame(corr_results)
    corr_df.to_csv(RESULTS_DIR / "proliferation_estrogen_correlation.csv", index=False)

    # Scatter plot: one panel per treatment
    ncols = len(available)
    nrows = 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 4), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).flatten()

    for i, treatment in enumerate(available):
        ax = axes[i]
        mask = adata.obs['treatment'] == treatment
        prolif = adata.obs.loc[mask, 'PROLIFERATION_score']
        er = adata.obs.loc[mask, er_score_col]

        ax.scatter(prolif, er, c=TREATMENT_COLORS.get(treatment, '#333'),
                  s=1, alpha=0.3, rasterized=True)

        # Add correlation text
        rho = corr_df.loc[corr_df['treatment'] == treatment, 'spearman_rho'].values[0]
        ax.text(0.05, 0.95, f'rho={rho:.3f}', transform=ax.transAxes,
               fontsize=10, va='top', fontweight='bold')
        ax.set_title(treatment, fontsize=10)

        if i % ncols == 0:
            ax.set_ylabel(er_score_col.replace('_score', '').replace('_', ' '))
        if i >= (nrows - 1) * ncols:
            ax.set_xlabel('Proliferation Score')

    plt.suptitle('Proliferation-Estrogen Response Coupling by Treatment', fontsize=12, y=1.01)
    save_fig(FIG_DIR / "panel7_prolif_er_coupling.png")
    plt.close()
    log_msg(f"  Saved: panel7_prolif_er_coupling.png/.pdf")
    log_msg(f"  Saved: {RESULTS_DIR / 'proliferation_estrogen_correlation.csv'}")

    return corr_df


# =============================================================================
# MAIN
# =============================================================================

def main():
    log_msg("=" * 60)
    log_msg("Proliferation & Quiescence Analysis")
    log_msg("=" * 60)

    # Load data
    log_msg(f"\nLoading data from {INPUT_H5AD}")
    check_file_exists(INPUT_H5AD, "Pathway-scored h5ad")
    adata = sc.read_h5ad(INPUT_H5AD)
    standardize_treatments(adata)
    log_msg(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    log_msg(f"Treatments: {adata.obs['treatment'].value_counts().to_dict()}")

    # Verify PROLIFERATION_score exists
    if 'PROLIFERATION_score' not in adata.obs.columns:
        raise ValueError("PROLIFERATION_score not found in adata.obs. Run script 07 first.")

    # Step 1: Score quiescence
    log_msg("\n" + "=" * 40)
    log_msg("STEP 1: Score quiescence")
    log_msg("=" * 40)
    adata = score_quiescence(adata)

    # Step 2: Classify cycling
    log_msg("\n" + "=" * 40)
    log_msg("STEP 2: Classify cycling vs non-cycling")
    log_msg("=" * 40)
    adata = classify_cycling(adata)

    # Step 3: Statistical comparisons
    log_msg("\n" + "=" * 40)
    log_msg("STEP 3: Statistical comparisons")
    log_msg("=" * 40)
    score_df, cycling_df = run_all_statistics(adata)

    # Save statistics
    score_df.to_csv(RESULTS_DIR / "score_comparisons.csv", index=False)
    cycling_df.to_csv(RESULTS_DIR / "cycling_fraction_comparisons.csv", index=False)
    log_msg(f"  Saved: {RESULTS_DIR / 'score_comparisons.csv'}")
    log_msg(f"  Saved: {RESULTS_DIR / 'cycling_fraction_comparisons.csv'}")

    # Step 4: Core visualizations (Panels 1-4)
    log_msg("\n" + "=" * 40)
    log_msg("STEP 4: Core visualizations")
    log_msg("=" * 40)
    plot_score_violins(adata)
    plot_cycling_fraction(adata)
    plot_key_markers(adata)
    plot_score_ridges(adata)

    # Step 5: Heterogeneity visualizations (Panels 5-7)
    log_msg("\n" + "=" * 40)
    log_msg("STEP 5: Heterogeneity visualizations")
    log_msg("=" * 40)
    plot_umap_overlays(adata)
    cluster_df = plot_cluster_cycling_heatmap(adata)
    corr_df = plot_proliferation_estrogen_coupling(adata)

    # Step 6: Save output
    log_msg("\n" + "=" * 40)
    log_msg("STEP 6: Save output")
    log_msg("=" * 40)
    adata.write_h5ad(OUTPUT_H5AD)
    log_msg(f"  Saved: {OUTPUT_H5AD}")

    # Summary
    log_msg("\n" + "=" * 60)
    log_msg("ANALYSIS COMPLETE")
    log_msg("=" * 60)
    log_msg(f"  Total cells: {adata.n_obs}")
    log_msg(f"  Cycling threshold: {adata.uns['cycling_threshold']:.3f} ({adata.uns['cycling_method']})")
    cycling_counts = adata.obs['cycling_status'].value_counts()
    log_msg(f"  Cycling: {cycling_counts.get('cycling', 0)} | Non-cycling: {cycling_counts.get('non-cycling', 0)}")
    log_msg(f"\n  Outputs:")
    log_msg(f"    h5ad: {OUTPUT_H5AD}")
    log_msg(f"    Results: {RESULTS_DIR}/")
    log_msg(f"    Figures: {FIG_DIR}/")

    return adata


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}")
        import traceback
        traceback.print_exc()
        raise
