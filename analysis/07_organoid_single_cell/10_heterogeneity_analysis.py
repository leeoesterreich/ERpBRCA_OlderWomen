#!/usr/bin/env python
# analysis/07_organoid_single_cell/10_heterogeneity_analysis.py
# Single-cell heterogeneity and biological axis analysis
#
# Explores within-condition heterogeneity, between-condition overlap in
# transcriptomic space, and identifies the core biological axis driving
# variation across treatment conditions.
#
# Analyses:
#   1. Per-condition UMAP density (KDE contours)
#   2. Neighborhood mixing analysis (entropy, mixing matrix)
#   3. Within-condition heterogeneity (dispersion, sub-clustering)
#   4. Cluster composition & enrichment (Fisher's exact)
#   5. Core biological axis (pathway PCA)
#   6. Estrogen response continuum (KDE overlap, composition)
#
# Inputs:
#   - data/processed/mechanism_explored.h5ad (from 08_inhibitor_mechanism_exploration.py)
#
# Outputs:
#   - data/processed/heterogeneity_analyzed.h5ad
#   - outputs/heterogeneity/mixing_matrix.csv
#   - outputs/heterogeneity/within_condition_heterogeneity.csv
#   - outputs/heterogeneity/cluster_treatment_contingency.csv
#   - outputs/heterogeneity/cluster_treatment_enrichment.csv
#   - outputs/heterogeneity/biological_axis_loadings.csv
#   - outputs/heterogeneity/estrogen_continuum_overlap.csv
#   - figures/heterogeneity/01_umap_density_per_condition.png/.pdf
#   - figures/heterogeneity/02_neighborhood_mixing.png/.pdf
#   - figures/heterogeneity/03_within_condition_dispersion.png/.pdf
#   - figures/heterogeneity/03b_umap_hulls_centroids.png/.pdf
#   - figures/heterogeneity/04_cluster_composition.png/.pdf
#   - figures/heterogeneity/05_biological_axis.png/.pdf
#   - figures/heterogeneity/06a_estrogen_continuum_density.png/.pdf
#   - figures/heterogeneity/06b_estrogen_continuum_composition.png/.pdf
#   - figures/heterogeneity/06c_estrogen_overlap_umap.png/.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _figure_config import save_fig
import seaborn as sns
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from scipy.spatial import ConvexHull
from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_samples
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
INPUT_H5AD = h5ad_path("mechanism_explored.h5ad")
OUTPUT_H5AD = h5ad_path("heterogeneity_analyzed.h5ad")
RESULTS_DIR = OUTPUT_DIR / "heterogeneity"
FIG_DIR = FIGURES_DIR / "heterogeneity"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Treatment ordering and colors
TREATMENT_ORDER = ['Vehicle', 'E1', 'E1+fulv', 'E1+HSD17B7i', 'E2', 'E2+fulv', 'E2+HSD17B7i']
TREATMENT_COLORS = {
    'Vehicle': '#999999',
    'E1': '#E41A1C',
    'E1+fulv': '#FF7F00',
    'E1+HSD17B7i': '#984EA3',
    'E2': '#377EB8',
    'E2+fulv': '#4DAF4A',
    'E2+HSD17B7i': '#A65628',
}

# Score columns for biological axis analysis
PATHWAY_SCORE_COLS = [
    'ESTROGEN_RESPONSE_EARLY_score', 'ESTROGEN_RESPONSE_LATE_score',
    'E2F_TARGETS_score', 'G2M_CHECKPOINT_score', 'MYC_TARGETS_score',
    'PROLIFERATION_score', 'STEROID_BIOSYNTHESIS_score', 'ER_TARGETS_DIRECT_score',
    'ER_PROLIFERATIVE_score', 'ER_TRANSCRIPTIONAL_score',
    'ER_METABOLIC_score', 'ER_SIGNALING_CROSSTALK_score',
    'S_score', 'G2M_score',
]


def load_data():
    """Load preprocessed h5ad."""
    print(f"[{datetime.now()}] Loading data from {INPUT_H5AD}")
    check_file_exists(INPUT_H5AD, "Mechanism-explored h5ad")
    adata = sc.read_h5ad(INPUT_H5AD)
    standardize_treatments(adata)
    print(f"  Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    # Ensure treatment is categorical with correct order
    adata.obs['treatment'] = pd.Categorical(
        adata.obs['treatment'], categories=TREATMENT_ORDER, ordered=True
    )
    return adata


# =============================================================================
# MODULE 1: Per-Condition UMAP Density
# =============================================================================
def module1_umap_density(adata):
    """KDE density contours per condition on UMAP."""
    print(f"\n[{datetime.now()}] === Module 1: Per-Condition UMAP Density ===")

    umap = adata.obsm['X_umap']
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    # Panel 0: all cells reference
    ax = axes[0]
    ax.scatter(umap[:, 0], umap[:, 1], s=0.5, alpha=0.1, c='grey', rasterized=True)
    ax.set_title('All Cells', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP1', fontsize=10)
    ax.set_ylabel('UMAP2', fontsize=10)

    for i, treat in enumerate(TREATMENT_ORDER):
        ax = axes[i + 1]
        mask = adata.obs['treatment'] == treat
        color = TREATMENT_COLORS[treat]

        # Gray background
        ax.scatter(umap[:, 0], umap[:, 1], s=0.3, alpha=0.05, c='grey', rasterized=True)

        # Condition points
        pts = umap[mask]
        ax.scatter(pts[:, 0], pts[:, 1], s=0.5, alpha=0.15, c=color, rasterized=True)

        # KDE contours
        try:
            # Subsample for KDE if too many cells
            if pts.shape[0] > 5000:
                idx = np.random.choice(pts.shape[0], 5000, replace=False)
                kde_pts = pts[idx]
            else:
                kde_pts = pts
            kde = stats.gaussian_kde(kde_pts.T, bw_method=0.3)
            xmin, xmax = umap[:, 0].min() - 1, umap[:, 0].max() + 1
            ymin, ymax = umap[:, 1].min() - 1, umap[:, 1].max() + 1
            xx, yy = np.meshgrid(np.linspace(xmin, xmax, 100), np.linspace(ymin, ymax, 100))
            zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
            ax.contour(xx, yy, zz, levels=5, colors=color, linewidths=1.0, alpha=0.8)
        except Exception as e:
            print(f"  KDE failed for {treat}: {e}")

        ax.set_title(f'{treat} (n={mask.sum():,})', fontsize=11, fontweight='bold')
        ax.set_xlabel('UMAP1', fontsize=10)
        ax.set_ylabel('UMAP2', fontsize=10)

    plt.suptitle('Per-Condition UMAP Density', fontsize=14, fontweight='bold', y=1.02)
    save_fig(FIG_DIR / "01_umap_density_per_condition.png")
    plt.close()
    print(f"  Saved UMAP density figure")


# =============================================================================
# MODULE 2: Neighborhood Mixing Analysis
# =============================================================================
def module2_neighborhood_mixing(adata):
    """Compute per-cell mixing entropy and condition mixing matrix."""
    print(f"\n[{datetime.now()}] === Module 2: Neighborhood Mixing Analysis ===")

    conn = adata.obsp['connectivities']
    treatments = adata.obs['treatment'].values
    treat_cats = TREATMENT_ORDER
    n_treats = len(treat_cats)
    treat_to_idx = {t: i for i, t in enumerate(treat_cats)}
    treat_idx = np.array([treat_to_idx[t] for t in treatments])

    n_cells = adata.n_obs
    mixing_fracs = np.zeros((n_cells, n_treats), dtype=np.float32)

    print(f"  Computing neighbor composition for {n_cells:,} cells...")
    # Process in chunks for memory efficiency
    chunk_size = 5000
    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = conn[start:end]
        if issparse(chunk):
            rows, cols = chunk.nonzero()
        else:
            rows, cols = np.where(chunk > 0)
        for r, c in zip(rows, cols):
            mixing_fracs[start + r, treat_idx[c]] += 1

    # Normalize rows
    row_sums = mixing_fracs.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    mixing_fracs /= row_sums

    # Shannon entropy per cell
    with np.errstate(divide='ignore', invalid='ignore'):
        log_fracs = np.where(mixing_fracs > 0, np.log2(mixing_fracs), 0)
        entropy = -np.sum(mixing_fracs * log_fracs, axis=1)

    adata.obs['mixing_entropy'] = entropy
    print(f"  Entropy range: [{entropy.min():.3f}, {entropy.max():.3f}]")

    # --- Mixing matrix ---
    mixing_matrix = pd.DataFrame(np.zeros((n_treats, n_treats)),
                                 index=treat_cats, columns=treat_cats)
    for i, treat in enumerate(treat_cats):
        mask = treat_idx == i
        mixing_matrix.iloc[i] = mixing_fracs[mask].mean(axis=0)

    mixing_matrix.to_csv(RESULTS_DIR / "mixing_matrix.csv")
    print(f"  Saved mixing matrix")

    # --- Figures ---
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))

    # UMAP colored by entropy
    ax = axes[0]
    sc_plot = ax.scatter(
        adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
        c=entropy, cmap='viridis', s=0.5, alpha=0.3, rasterized=True
    )
    plt.colorbar(sc_plot, ax=ax, label='Mixing Entropy')
    ax.set_title('Neighborhood Mixing Entropy', fontsize=14, fontweight='bold')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    # Mixing matrix heatmap
    ax = axes[1]
    sns.heatmap(mixing_matrix, annot=True, fmt='.2f', cmap='YlOrRd',
                ax=ax, square=True, cbar_kws={'label': 'Neighbor Fraction'})
    ax.set_title('Condition Mixing Matrix', fontsize=14, fontweight='bold')
    ax.set_xlabel('Neighbor Condition')
    ax.set_ylabel('Cell Condition')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=11)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=45, ha='right', fontsize=11)

    # Entropy violin by condition
    ax = axes[2]
    plot_df = pd.DataFrame({
        'Mixing Entropy': entropy,
        'Treatment': treatments
    })
    sns.violinplot(data=plot_df, x='Treatment', y='Mixing Entropy',
                   order=TREATMENT_ORDER, palette=TREATMENT_COLORS,
                   ax=ax, inner='box', cut=0, scale='width')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_title('Mixing Entropy by Condition', fontsize=14, fontweight='bold')

    save_fig(FIG_DIR / "02_neighborhood_mixing.png")
    plt.close()
    print(f"  Saved mixing analysis figure")


# =============================================================================
# MODULE 3: Within-Condition Heterogeneity
# =============================================================================
def module3_within_condition_heterogeneity(adata):
    """Measure dispersion, silhouette, and sub-clustering per condition."""
    print(f"\n[{datetime.now()}] === Module 3: Within-Condition Heterogeneity ===")

    umap = adata.obsm['X_umap']
    pca = adata.obsm['X_pca'][:, :30]
    treatments = adata.obs['treatment'].values

    results = []
    sub_cluster_all = pd.Series('', index=adata.obs.index, dtype=str)

    # Silhouette scores in PCA space (subsample for speed)
    print(f"  Computing silhouette scores...")
    n_sil = min(10000, adata.n_obs)
    sil_idx = np.random.choice(adata.n_obs, n_sil, replace=False)
    sil_scores = silhouette_samples(pca[sil_idx], treatments[sil_idx])
    sil_df = pd.DataFrame({'treatment': treatments[sil_idx], 'silhouette': sil_scores})
    sil_per_treat = sil_df.groupby('treatment')['silhouette'].mean()

    for treat in TREATMENT_ORDER:
        mask = treatments == treat
        pts_umap = umap[mask]
        pts_pca = pca[mask]
        n_cells = mask.sum()

        # UMAP dispersion: mean distance to centroid
        centroid_umap = pts_umap.mean(axis=0)
        dist_to_centroid = np.sqrt(((pts_umap - centroid_umap) ** 2).sum(axis=1))
        umap_dispersion = dist_to_centroid.mean()

        # Convex hull area
        try:
            hull = ConvexHull(pts_umap)
            hull_area = hull.volume  # 2D: volume = area
        except Exception:
            hull_area = np.nan

        # PCA dispersion: mean pairwise distance (subsampled)
        n_sub = min(1000, n_cells)
        sub_idx = np.random.choice(n_cells, n_sub, replace=False)
        from scipy.spatial.distance import pdist
        pca_dists = pdist(pts_pca[sub_idx])
        pca_dispersion = pca_dists.mean()

        # Silhouette
        sil = sil_per_treat.get(treat, np.nan)

        # Sub-clustering within condition
        adata_sub = adata[mask].copy()
        sc.pp.neighbors(adata_sub, n_pcs=30, use_rep='X_pca')
        sc.tl.leiden(adata_sub, resolution=0.3, key_added='sub_leiden')
        n_subclusters = adata_sub.obs['sub_leiden'].nunique()
        sub_cluster_all[mask] = treat + '_' + adata_sub.obs['sub_leiden'].astype(str).values

        results.append({
            'treatment': treat,
            'n_cells': n_cells,
            'umap_dispersion': umap_dispersion,
            'hull_area': hull_area,
            'pca_dispersion': pca_dispersion,
            'silhouette': sil,
            'n_subclusters': n_subclusters,
        })
        print(f"  {treat}: dispersion={umap_dispersion:.2f}, hull={hull_area:.1f}, "
              f"sil={sil:.3f}, subclusters={n_subclusters}")

    adata.obs['sub_cluster'] = sub_cluster_all
    results_df = pd.DataFrame(results)
    results_df.to_csv(RESULTS_DIR / "within_condition_heterogeneity.csv", index=False)

    # --- Figures ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Dispersion bar chart
    ax = axes[0]
    x = np.arange(len(TREATMENT_ORDER))
    colors = [TREATMENT_COLORS[t] for t in TREATMENT_ORDER]
    ax.bar(x, results_df['umap_dispersion'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(TREATMENT_ORDER, rotation=45, ha='right')
    ax.set_ylabel('Mean Dist. to Centroid (UMAP)')
    ax.set_title('UMAP Dispersion', fontsize=11, fontweight='bold')

    # PCA dispersion
    ax = axes[1]
    ax.bar(x, results_df['pca_dispersion'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(TREATMENT_ORDER, rotation=45, ha='right')
    ax.set_ylabel('Mean Pairwise Distance (PCA)')
    ax.set_title('PCA-Space Dispersion', fontsize=11, fontweight='bold')

    # Silhouette scores
    ax = axes[2]
    ax.bar(x, results_df['silhouette'], color=colors, edgecolor='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(TREATMENT_ORDER, rotation=45, ha='right')
    ax.set_ylabel('Mean Silhouette Score')
    ax.set_title('Condition Separation (Silhouette)', fontsize=11, fontweight='bold')
    ax.axhline(0, color='black', linewidth=0.5, linestyle='--')

    save_fig(FIG_DIR / "03_within_condition_dispersion.png")
    plt.close()

    # UMAP with convex hulls and centroids
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(umap[:, 0], umap[:, 1], s=0.3, alpha=0.05, c='grey', rasterized=True)
    for treat in TREATMENT_ORDER:
        mask = treatments == treat
        pts = umap[mask]
        color = TREATMENT_COLORS[treat]
        centroid = pts.mean(axis=0)
        ax.scatter(*centroid, s=100, c=color, edgecolors='black', linewidths=1.5,
                   zorder=5, marker='*', label=treat)
        try:
            hull = ConvexHull(pts)
            for simplex in hull.simplices:
                ax.plot(pts[simplex, 0], pts[simplex, 1], c=color, alpha=0.5, linewidth=1)
        except Exception:
            pass
    ax.legend(markerscale=1.5, fontsize=10, loc='best')
    ax.set_title('Condition Centroids & Convex Hulls', fontsize=12, fontweight='bold')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    save_fig(FIG_DIR / "03b_umap_hulls_centroids.png")
    plt.close()
    print(f"  Saved heterogeneity figures")


# =============================================================================
# MODULE 4: Cluster Composition & Enrichment
# =============================================================================
def module4_cluster_composition(adata):
    """Contingency analysis of Leiden clusters x treatment with Fisher's exact."""
    print(f"\n[{datetime.now()}] === Module 4: Cluster Composition & Enrichment ===")

    clusters = adata.obs['leiden'].values
    treatments = adata.obs['treatment'].values
    cluster_cats = sorted(adata.obs['leiden'].cat.categories.tolist(),
                          key=lambda x: int(x) if x.isdigit() else x)

    # Contingency table
    contingency = pd.crosstab(
        pd.Categorical(clusters, categories=cluster_cats),
        pd.Categorical(treatments, categories=TREATMENT_ORDER)
    )
    contingency.to_csv(RESULTS_DIR / "cluster_treatment_contingency.csv")

    # Fisher's exact test per (cluster, condition) pair
    enrichment_results = []
    total = len(clusters)
    for cl in cluster_cats:
        in_cluster = clusters == cl
        for treat in TREATMENT_ORDER:
            in_treat = treatments == treat
            a = np.sum(in_cluster & in_treat)
            b = np.sum(in_cluster & ~in_treat)
            c = np.sum(~in_cluster & in_treat)
            d = np.sum(~in_cluster & ~in_treat)
            odds, pval = stats.fisher_exact([[a, b], [c, d]])
            enrichment_results.append({
                'cluster': cl, 'treatment': treat,
                'observed': a, 'odds_ratio': odds, 'pvalue': pval,
                'pct_of_cluster': a / (a + b) * 100 if (a + b) > 0 else 0,
            })

    enrich_df = pd.DataFrame(enrichment_results)
    # BH correction
    from statsmodels.stats.multitest import multipletests
    _, enrich_df['padj'], _, _ = multipletests(enrich_df['pvalue'], method='fdr_bh')
    enrich_df['sig'] = enrich_df['padj'] < 0.05
    enrich_df.to_csv(RESULTS_DIR / "cluster_treatment_enrichment.csv", index=False)
    print(f"  {enrich_df['sig'].sum()} significant enrichments (FDR < 0.05)")

    # --- Figures ---
    fig, axes = plt.subplots(1, 3, figsize=(22, 6))

    # Stacked bar chart (proportion)
    ax = axes[0]
    prop = contingency.div(contingency.sum(axis=1), axis=0)
    prop.plot(kind='bar', stacked=True, ax=ax,
              color=[TREATMENT_COLORS[t] for t in TREATMENT_ORDER],
              edgecolor='white', linewidth=0.3)
    ax.set_ylabel('Proportion')
    ax.set_xlabel('Leiden Cluster')
    ax.set_title('Cluster Composition by Treatment', fontsize=11, fontweight='bold')
    ax.legend(title='Treatment', bbox_to_anchor=(1.0, 1.0), fontsize=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    # Enrichment heatmap (log2 odds ratio, significance marked)
    ax = axes[1]
    pivot_or = enrich_df.pivot(index='cluster', columns='treatment', values='odds_ratio')
    pivot_or = pivot_or.reindex(index=cluster_cats, columns=TREATMENT_ORDER)
    log2_or = np.log2(pivot_or.clip(lower=0.01))

    pivot_sig = enrich_df.pivot(index='cluster', columns='treatment', values='sig')
    pivot_sig = pivot_sig.reindex(index=cluster_cats, columns=TREATMENT_ORDER)

    sns.heatmap(log2_or, cmap='RdBu_r', center=0, ax=ax,
                annot=True, fmt='.1f', annot_kws={'fontsize': 10},
                cbar_kws={'label': 'log2(Odds Ratio)'})
    # Mark significant cells
    for i in range(len(cluster_cats)):
        for j in range(len(TREATMENT_ORDER)):
            if pivot_sig.iloc[i, j]:
                ax.text(j + 0.5, i + 0.85, '*', ha='center', va='center',
                        fontsize=12, fontweight='bold', color='black')
    ax.set_title('Enrichment (log2 OR, * = FDR<0.05)', fontsize=11, fontweight='bold')
    ax.set_xlabel('Treatment')
    ax.set_ylabel('Cluster')

    # Annotated UMAP
    ax = axes[2]
    umap = adata.obsm['X_umap']
    # Color by enrichment: which treatment is most enriched in each cluster
    dominant_treat = enrich_df.loc[enrich_df.groupby('cluster')['odds_ratio'].idxmax()]
    cluster_color_map = {row['cluster']: TREATMENT_COLORS[row['treatment']]
                         for _, row in dominant_treat.iterrows()}
    cell_colors = [cluster_color_map.get(c, '#cccccc') for c in clusters]
    ax.scatter(umap[:, 0], umap[:, 1], c=cell_colors, s=0.5, alpha=0.2, rasterized=True)
    # Add cluster labels
    for cl in cluster_cats:
        mask = clusters == cl
        cx, cy = umap[mask, 0].mean(), umap[mask, 1].mean()
        ax.text(cx, cy, cl, fontsize=10, fontweight='bold',
                ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7))
    ax.set_title('Clusters Colored by Dominant Treatment', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    save_fig(FIG_DIR / "04_cluster_composition.png")
    plt.close()
    print(f"  Saved cluster composition figures")


# =============================================================================
# MODULE 5: Core Biological Axis (Pathway PCA)
# =============================================================================
def module5_biological_axis(adata):
    """PCA on pathway scores to find dominant biological axis."""
    print(f"\n[{datetime.now()}] === Module 5: Core Biological Axis ===")

    # Gather score columns that exist
    available_cols = [c for c in PATHWAY_SCORE_COLS if c in adata.obs.columns]
    print(f"  Using {len(available_cols)} score columns")

    score_df = adata.obs[available_cols].copy()
    # Standardize
    scaler = StandardScaler()
    scaled = scaler.fit_transform(score_df.values)

    # PCA
    pca = PCA(n_components=min(5, len(available_cols)))
    pcs = pca.fit_transform(scaled)

    adata.obs['biological_axis_PC1'] = pcs[:, 0]
    adata.obs['biological_axis_PC2'] = pcs[:, 1]

    # Loadings
    loadings = pd.DataFrame(
        pca.components_.T,
        index=available_cols,
        columns=[f'PC{i+1}' for i in range(pca.n_components_)]
    )
    loadings.to_csv(RESULTS_DIR / "biological_axis_loadings.csv")

    var_explained = pca.explained_variance_ratio_
    print(f"  PC1 explains {var_explained[0]*100:.1f}% variance")
    print(f"  Top PC1 loadings:")
    sorted_load = loadings['PC1'].abs().sort_values(ascending=False)
    for gene, val in sorted_load.head(5).items():
        sign = '+' if loadings.loc[gene, 'PC1'] > 0 else '-'
        print(f"    {sign}{gene}: {loadings.loc[gene, 'PC1']:.3f}")

    # --- Figures ---
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    # PC1 loadings bar chart
    ax = axes[0, 0]
    load_sorted = loadings['PC1'].sort_values()
    colors_load = ['#E41A1C' if v > 0 else '#377EB8' for v in load_sorted.values]
    # Clean up names for display
    labels = [c.replace('_score', '').replace('_', ' ') for c in load_sorted.index]
    ax.barh(range(len(load_sorted)), load_sorted.values, color=colors_load, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('PC1 Loading')
    ax.set_title(f'PC1 Loadings ({var_explained[0]*100:.1f}% var)', fontsize=11, fontweight='bold')
    ax.axvline(0, color='black', linewidth=0.5)

    # Scree plot
    ax = axes[0, 1]
    ax.bar(range(1, len(var_explained) + 1), var_explained * 100,
           color='steelblue', edgecolor='black', linewidth=0.5)
    ax.plot(range(1, len(var_explained) + 1), np.cumsum(var_explained) * 100,
            'ro-', markersize=6)
    ax.set_xlabel('PC')
    ax.set_ylabel('Variance Explained (%)')
    ax.set_title('Scree Plot', fontsize=11, fontweight='bold')

    # UMAP colored by PC1
    ax = axes[0, 2]
    vmin, vmax = np.percentile(pcs[:, 0], [2, 98])
    sc_plot = ax.scatter(
        adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
        c=pcs[:, 0], cmap='RdBu_r', s=0.5, alpha=0.3, vmin=vmin, vmax=vmax,
        rasterized=True
    )
    plt.colorbar(sc_plot, ax=ax, label='PC1 Score')
    ax.set_title('Biological Axis PC1', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    # Violin of PC1 by condition
    ax = axes[1, 0]
    plot_df = pd.DataFrame({
        'PC1': pcs[:, 0],
        'Treatment': adata.obs['treatment'].values
    })
    sns.violinplot(data=plot_df, x='Treatment', y='PC1',
                   order=TREATMENT_ORDER, palette=TREATMENT_COLORS,
                   ax=ax, inner='box', cut=0, scale='width')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_title('PC1 by Condition', fontsize=11, fontweight='bold')
    ax.set_ylabel('Biological Axis PC1')

    # PC1 vs PC2 scatter colored by condition
    ax = axes[1, 1]
    for treat in TREATMENT_ORDER:
        mask = adata.obs['treatment'].values == treat
        ax.scatter(pcs[mask, 0], pcs[mask, 1], s=1, alpha=0.15,
                   c=TREATMENT_COLORS[treat], label=treat, rasterized=True)
    ax.set_xlabel(f'PC1 ({var_explained[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({var_explained[1]*100:.1f}%)')
    ax.set_title('Pathway PCA Space', fontsize=11, fontweight='bold')
    ax.legend(markerscale=5, fontsize=10, loc='best')

    # PC2 loadings
    ax = axes[1, 2]
    load2_sorted = loadings['PC2'].sort_values()
    colors2 = ['#E41A1C' if v > 0 else '#377EB8' for v in load2_sorted.values]
    labels2 = [c.replace('_score', '').replace('_', ' ') for c in load2_sorted.index]
    ax.barh(range(len(load2_sorted)), load2_sorted.values, color=colors2, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(labels2)))
    ax.set_yticklabels(labels2, fontsize=10)
    ax.set_xlabel('PC2 Loading')
    ax.set_title(f'PC2 Loadings ({var_explained[1]*100:.1f}% var)', fontsize=11, fontweight='bold')
    ax.axvline(0, color='black', linewidth=0.5)

    plt.suptitle('Core Biological Axis Analysis', fontsize=14, fontweight='bold', y=1.02)
    save_fig(FIG_DIR / "05_biological_axis.png")
    plt.close()
    print(f"  Saved biological axis figures")


# =============================================================================
# MODULE 6: Estrogen Response Continuum
# =============================================================================
def module6_estrogen_continuum(adata):
    """Analyze the estrogen response as a continuous spectrum."""
    print(f"\n[{datetime.now()}] === Module 6: Estrogen Response Continuum ===")

    # Composite estrogen score
    er_early = adata.obs['ESTROGEN_RESPONSE_EARLY_score'].values
    er_late = adata.obs['ESTROGEN_RESPONSE_LATE_score'].values
    continuum = (er_early + er_late) / 2
    adata.obs['estrogen_continuum'] = continuum
    print(f"  Continuum range: [{continuum.min():.3f}, {continuum.max():.3f}]")

    treatments = adata.obs['treatment'].values

    # --- Figure 1: KDE density curves ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))

    e1_group = ['Vehicle', 'E1', 'E1+fulv', 'E1+HSD17B7i']
    e2_group = ['Vehicle', 'E2', 'E2+fulv', 'E2+HSD17B7i']

    for ax, group, title in [(axes[0], e1_group, 'E1 Group'),
                              (axes[1], e2_group, 'E2 Group')]:
        x_range = np.linspace(continuum.min(), continuum.max(), 200)
        for treat in group:
            mask = treatments == treat
            kde = stats.gaussian_kde(continuum[mask], bw_method=0.2)
            ax.plot(x_range, kde(x_range), label=treat, color=TREATMENT_COLORS[treat],
                    linewidth=2)
            ax.fill_between(x_range, kde(x_range), alpha=0.1, color=TREATMENT_COLORS[treat])
        ax.set_xlabel('Estrogen Response Continuum')
        ax.set_ylabel('Density')
        ax.set_title(f'{title}: Estrogen Response Distribution', fontsize=11, fontweight='bold')
        ax.legend(fontsize=10)

    save_fig(FIG_DIR / "06a_estrogen_continuum_density.png")
    plt.close()

    # --- Figure 2: Sliding window composition ---
    fig, ax = plt.subplots(figsize=(12, 5))
    n_windows = 50
    window_centers = np.linspace(np.percentile(continuum, 2), np.percentile(continuum, 98), n_windows)
    window_width = (window_centers[1] - window_centers[0]) * 2
    composition = np.zeros((n_windows, len(TREATMENT_ORDER)))
    for i, center in enumerate(window_centers):
        in_window = np.abs(continuum - center) < window_width / 2
        for j, treat in enumerate(TREATMENT_ORDER):
            composition[i, j] = np.sum(in_window & (treatments == treat))
    # Normalize
    row_sums = composition.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    composition /= row_sums

    ax.stackplot(window_centers, composition.T,
                 labels=TREATMENT_ORDER,
                 colors=[TREATMENT_COLORS[t] for t in TREATMENT_ORDER],
                 alpha=0.8)
    ax.set_xlabel('Estrogen Response Continuum')
    ax.set_ylabel('Proportion')
    ax.set_title('Condition Composition Along Estrogen Continuum', fontsize=12, fontweight='bold')
    ax.legend(loc='upper left', fontsize=10)
    ax.set_xlim(window_centers[0], window_centers[-1])
    ax.set_ylim(0, 1)
    save_fig(FIG_DIR / "06b_estrogen_continuum_composition.png")
    plt.close()

    # --- Figure 3: Pairwise KDE overlap ---
    # Compute overlap as integral of min(kde1, kde2)
    x_eval = np.linspace(continuum.min() - 0.5, continuum.max() + 0.5, 500)
    kdes = {}
    for treat in TREATMENT_ORDER:
        mask = treatments == treat
        kdes[treat] = stats.gaussian_kde(continuum[mask], bw_method=0.2)(x_eval)

    overlap_matrix = pd.DataFrame(np.zeros((len(TREATMENT_ORDER), len(TREATMENT_ORDER))),
                                  index=TREATMENT_ORDER, columns=TREATMENT_ORDER)
    dx = x_eval[1] - x_eval[0]
    for i, t1 in enumerate(TREATMENT_ORDER):
        for j, t2 in enumerate(TREATMENT_ORDER):
            overlap = np.sum(np.minimum(kdes[t1], kdes[t2])) * dx
            overlap_matrix.iloc[i, j] = overlap

    overlap_matrix.to_csv(RESULTS_DIR / "estrogen_continuum_overlap.csv")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Overlap heatmap
    ax = axes[0]
    sns.heatmap(overlap_matrix, annot=True, fmt='.2f', cmap='YlOrRd',
                ax=ax, square=True, vmin=0, vmax=1,
                cbar_kws={'label': 'KDE Overlap'})
    ax.set_title('Pairwise Estrogen Response Overlap', fontsize=11, fontweight='bold')

    # UMAP colored by continuum
    ax = axes[1]
    vmin, vmax = np.percentile(continuum, [2, 98])
    sc_plot = ax.scatter(
        adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
        c=continuum, cmap='RdYlBu_r', s=0.5, alpha=0.3, vmin=vmin, vmax=vmax,
        rasterized=True
    )
    plt.colorbar(sc_plot, ax=ax, label='Estrogen Response')
    ax.set_title('Estrogen Response Continuum', fontsize=11, fontweight='bold')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    save_fig(FIG_DIR / "06c_estrogen_overlap_umap.png")
    plt.close()
    print(f"  Saved estrogen continuum figures")


# =============================================================================
# MAIN
# =============================================================================
def main():
    print("=" * 70)
    print("SINGLE-CELL HETEROGENEITY & BIOLOGICAL AXIS ANALYSIS")
    print("=" * 70)

    np.random.seed(42)
    adata = load_data()

    module1_umap_density(adata)
    module2_neighborhood_mixing(adata)
    module3_within_condition_heterogeneity(adata)
    module4_cluster_composition(adata)
    module5_biological_axis(adata)
    module6_estrogen_continuum(adata)

    # Save annotated object
    print(f"\n[{datetime.now()}] Saving annotated AnnData to {OUTPUT_H5AD}")
    adata.write_h5ad(OUTPUT_H5AD)

    print(f"\n[{datetime.now()}] === COMPLETE ===")
    print(f"  Figures: {FIG_DIR}/")
    print(f"  Results: {RESULTS_DIR}/")
    print(f"  AnnData: {OUTPUT_H5AD}")


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}")
        import traceback
        traceback.print_exc()
        raise
