#!/usr/bin/env python
# analysis/07_organoid_single_cell/03_preprocess.py
# Preprocessing pipeline: normalization, HVG selection, PCA, UMAP, clustering
#
# Inputs:
#   - data/processed/qc_filtered.h5ad (from 02_qc.py)
#
# Outputs:
#   - data/processed/preprocessed.h5ad
#   - figures/umap/umap_by_sample.png/.pdf
#   - figures/umap/umap_by_treatment.png/.pdf
#   - figures/umap/umap_by_cluster.png/.pdf
#   - figures/umap/umap_by_phase.png/.pdf
#   - figures/umap/umap_epithelial_markers.png/.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scanpy as sc
import pandas as pd
import numpy as np
import os
from datetime import datetime

from _config import (
    FIGURES_DIR, h5ad_path, check_file_exists
)

# =============================================================================
# Configuration
# =============================================================================
# Input/Output paths
QC_FILTERED_FILE = h5ad_path("qc_filtered.h5ad")
PREPROCESSED_FILE = h5ad_path("preprocessed.h5ad")

# Figure directories
UMAP_FIGURES_DIR = FIGURES_DIR / "umap"

# Cell cycle gene sets (Tirosh 2016)
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN",
    "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
]

# Epithelial markers to plot
EPITHELIAL_MARKERS = [
    "KRT8", "KRT18", "EPCAM", "ESR1", "PGR", "KRT5", "KRT14", "TP63", "MKI67"
]

# Scanpy settings
sc.settings.verbosity = 2
sc.settings.figdir = str(UMAP_FIGURES_DIR)
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)


def log_message(msg):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}")


def load_qc_filtered_data():
    """Load QC filtered data from h5ad file.

    Returns
    -------
    adata : AnnData
        QC filtered AnnData object
    """
    log_message(f"Loading QC filtered data from {QC_FILTERED_FILE}")
    check_file_exists(QC_FILTERED_FILE, "QC filtered file")

    adata = sc.read_h5ad(QC_FILTERED_FILE)
    log_message(f"  Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    return adata


def store_raw_counts(adata):
    """Store raw counts in a layer before normalization.

    Parameters
    ----------
    adata : AnnData
        AnnData object with raw counts in X

    Returns
    -------
    adata : AnnData
        AnnData with counts stored in layers["counts"]
    """
    log_message("Storing raw counts in layer...")
    adata.layers["counts"] = adata.X.copy()
    log_message("  Stored raw counts in adata.layers['counts']")
    return adata


def filter_genes(adata, min_cells=3):
    """Filter genes expressed in fewer than min_cells.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    min_cells : int
        Minimum number of cells expressing a gene

    Returns
    -------
    adata : AnnData
        Filtered AnnData object
    """
    log_message(f"Filtering genes expressed in < {min_cells} cells...")
    n_genes_before = adata.n_vars

    sc.pp.filter_genes(adata, min_cells=min_cells)

    n_genes_after = adata.n_vars
    log_message(f"  Filtered from {n_genes_before} to {n_genes_after} genes "
                f"(removed {n_genes_before - n_genes_after})")

    return adata


def normalize_and_log(adata, target_sum=10000):
    """Normalize counts to target_sum per cell and log1p transform.

    Parameters
    ----------
    adata : AnnData
        AnnData object with counts
    target_sum : int
        Target counts per cell after normalization

    Returns
    -------
    adata : AnnData
        Normalized and log-transformed AnnData
    """
    log_message(f"Normalizing to {target_sum} counts per cell...")
    sc.pp.normalize_total(adata, target_sum=target_sum)

    log_message("Applying log1p transformation...")
    sc.pp.log1p(adata)

    log_message("  Normalization complete")
    return adata


def select_highly_variable_genes(adata, n_top_genes=3000):
    """Select highly variable genes using seurat_v3 method.

    Parameters
    ----------
    adata : AnnData
        AnnData object (normalized)
    n_top_genes : int
        Number of top HVGs to select

    Returns
    -------
    adata : AnnData
        AnnData with highly_variable column in var
    """
    log_message(f"Selecting {n_top_genes} highly variable genes (seurat_v3)...")

    # seurat_v3 uses counts layer
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor='seurat_v3',
        layer='counts',
        subset=False  # Keep all genes, just mark HVGs
    )

    n_hvg = adata.var['highly_variable'].sum()
    log_message(f"  Selected {n_hvg} highly variable genes")

    return adata


def scale_data(adata, max_value=10):
    """Scale data to unit variance and zero mean, clipping to max_value.

    Parameters
    ----------
    adata : AnnData
        Normalized AnnData object
    max_value : float
        Maximum value after scaling

    Returns
    -------
    adata : AnnData
        Scaled AnnData object
    """
    log_message(f"Scaling data (max_value={max_value})...")
    sc.pp.scale(adata, max_value=max_value)
    log_message("  Scaling complete")
    return adata


def run_pca(adata, n_comps=50):
    """Run PCA on highly variable genes.

    Parameters
    ----------
    adata : AnnData
        Scaled AnnData object
    n_comps : int
        Number of PCA components

    Returns
    -------
    adata : AnnData
        AnnData with PCA results in obsm['X_pca']
    """
    log_message(f"Running PCA ({n_comps} components) on HVGs...")

    # Run PCA using only HVGs
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)

    # Log variance explained
    var_explained = adata.uns['pca']['variance_ratio']
    cumsum_var = np.cumsum(var_explained)
    log_message(f"  Variance explained by top 10 PCs: {100 * cumsum_var[9]:.1f}%")
    log_message(f"  Variance explained by top 30 PCs: {100 * cumsum_var[29]:.1f}%")
    log_message(f"  Variance explained by all 50 PCs: {100 * cumsum_var[49]:.1f}%")

    return adata


def score_cell_cycle(adata, s_genes, g2m_genes):
    """Score cell cycle phases using Tirosh gene sets.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    s_genes : list
        List of S phase genes
    g2m_genes : list
        List of G2M phase genes

    Returns
    -------
    adata : AnnData
        AnnData with S_score, G2M_score, and phase in obs
    """
    log_message("Scoring cell cycle phases...")

    # Check which genes are present in the dataset
    all_genes = set(adata.var_names)

    s_genes_present = [g for g in s_genes if g in all_genes]
    g2m_genes_present = [g for g in g2m_genes if g in all_genes]

    s_genes_missing = set(s_genes) - set(s_genes_present)
    g2m_genes_missing = set(g2m_genes) - set(g2m_genes_present)

    log_message(f"  S genes: {len(s_genes_present)}/{len(s_genes)} present")
    if s_genes_missing:
        log_message(f"    Missing: {', '.join(sorted(s_genes_missing))}")

    log_message(f"  G2M genes: {len(g2m_genes_present)}/{len(g2m_genes)} present")
    if g2m_genes_missing:
        log_message(f"    Missing: {', '.join(sorted(g2m_genes_missing))}")

    # Score cell cycle
    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=s_genes_present,
        g2m_genes=g2m_genes_present
    )

    # Log phase distribution
    phase_counts = adata.obs['phase'].value_counts()
    log_message("  Phase distribution:")
    for phase, count in phase_counts.items():
        pct = 100 * count / adata.n_obs
        log_message(f"    {phase}: {count} cells ({pct:.1f}%)")

    return adata


def compute_neighbors(adata, n_neighbors=15, n_pcs=30):
    """Compute nearest neighbors for clustering and UMAP.

    Parameters
    ----------
    adata : AnnData
        AnnData with PCA
    n_neighbors : int
        Number of nearest neighbors
    n_pcs : int
        Number of PCs to use

    Returns
    -------
    adata : AnnData
        AnnData with neighbor graph in uns
    """
    log_message(f"Computing neighbors (n_neighbors={n_neighbors}, n_pcs={n_pcs})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    log_message("  Neighbor graph computed")
    return adata


def run_umap(adata):
    """Run UMAP dimensionality reduction.

    Parameters
    ----------
    adata : AnnData
        AnnData with neighbor graph

    Returns
    -------
    adata : AnnData
        AnnData with UMAP in obsm['X_umap']
    """
    log_message("Running UMAP...")
    sc.tl.umap(adata)
    log_message("  UMAP complete")
    return adata


def run_leiden_clustering(adata, resolution=0.5):
    """Run Leiden clustering.

    Parameters
    ----------
    adata : AnnData
        AnnData with neighbor graph
    resolution : float
        Resolution parameter for Leiden

    Returns
    -------
    adata : AnnData
        AnnData with 'leiden' cluster assignments in obs
    """
    log_message(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution)

    n_clusters = adata.obs['leiden'].nunique()
    log_message(f"  Found {n_clusters} clusters")

    # Log cluster sizes
    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    log_message("  Cluster sizes:")
    for cluster, count in cluster_counts.items():
        pct = 100 * count / adata.n_obs
        log_message(f"    Cluster {cluster}: {count} cells ({pct:.1f}%)")

    return adata


def plot_umap_by_category(adata, color_by, filename, title=None):
    """Plot UMAP colored by a categorical variable.

    Parameters
    ----------
    adata : AnnData
        AnnData with UMAP
    color_by : str
        Column in obs to color by
    filename : str
        Output filename (without extension)
    title : str, optional
        Plot title
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    sc.pl.umap(
        adata,
        color=color_by,
        ax=ax,
        show=False,
        title=title or color_by,
        legend_loc='right margin'
    )

    plt.tight_layout()
    out_path_png = UMAP_FIGURES_DIR / f"{filename}.png"
    out_path_pdf = UMAP_FIGURES_DIR / f"{filename}.pdf"
    plt.savefig(out_path_png, dpi=150, bbox_inches='tight')
    plt.savefig(out_path_pdf, bbox_inches='tight')
    plt.close()
    log_message(f"  Saved: {out_path_png}")


def save_umap_plots(adata):
    """Save UMAP plots for various categories.

    Parameters
    ----------
    adata : AnnData
        AnnData with UMAP and clustering
    """
    log_message("Saving UMAP plots...")

    # UMAP by sample
    plot_umap_by_category(
        adata,
        color_by='sample_id',
        filename='umap_by_sample',
        title='UMAP by Sample'
    )

    # UMAP by treatment
    plot_umap_by_category(
        adata,
        color_by='treatment',
        filename='umap_by_treatment',
        title='UMAP by Treatment'
    )

    # UMAP by cluster
    plot_umap_by_category(
        adata,
        color_by='leiden',
        filename='umap_by_cluster',
        title='UMAP by Leiden Cluster'
    )

    # UMAP by cell cycle phase
    plot_umap_by_category(
        adata,
        color_by='phase',
        filename='umap_by_phase',
        title='UMAP by Cell Cycle Phase'
    )


def plot_epithelial_markers(adata, markers):
    """Plot expression of epithelial markers on UMAP.

    Parameters
    ----------
    adata : AnnData
        AnnData with UMAP
    markers : list
        List of marker genes to plot
    """
    log_message("Plotting epithelial markers...")

    # Check which markers are present
    all_genes = set(adata.var_names)
    markers_present = [m for m in markers if m in all_genes]
    markers_missing = [m for m in markers if m not in all_genes]

    if markers_missing:
        log_message(f"  Warning: Missing markers: {', '.join(markers_missing)}")

    if not markers_present:
        log_message("  No markers found in dataset, skipping marker plots")
        return

    log_message(f"  Plotting {len(markers_present)} markers: {', '.join(markers_present)}")

    # Calculate grid dimensions
    n_markers = len(markers_present)
    ncols = min(3, n_markers)
    nrows = (n_markers + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))

    # Handle single row/column cases
    if n_markers == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Plot each marker
    for i, marker in enumerate(markers_present):
        sc.pl.umap(
            adata,
            color=marker,
            ax=axes[i],
            show=False,
            title=marker,
            color_map='viridis'
        )

    # Remove empty axes
    for j in range(len(markers_present), len(axes)):
        axes[j].axis('off')

    plt.tight_layout()
    out_path_png = UMAP_FIGURES_DIR / "umap_epithelial_markers.png"
    out_path_pdf = UMAP_FIGURES_DIR / "umap_epithelial_markers.pdf"
    plt.savefig(out_path_png, dpi=150, bbox_inches='tight')
    plt.savefig(out_path_pdf, bbox_inches='tight')
    plt.close()
    log_message(f"  Saved: {out_path_png}")


def save_preprocessed_data(adata):
    """Save preprocessed AnnData to h5ad file.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData object
    """
    log_message(f"Saving preprocessed data to {PREPROCESSED_FILE}")
    adata.write_h5ad(PREPROCESSED_FILE)
    log_message(f"  Saved: {adata.n_obs} cells, {adata.n_vars} genes")


def main():
    """Main preprocessing pipeline."""
    log_message("=" * 60)
    log_message("Starting preprocessing pipeline for organoid scRNA-seq")
    log_message("=" * 60)

    # Ensure output directories exist
    UMAP_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Load QC filtered data
    log_message("\n" + "=" * 40)
    log_message("STEP 1: Load QC filtered data")
    log_message("=" * 40)
    adata = load_qc_filtered_data()

    # Step 2: Store raw counts in layer
    log_message("\n" + "=" * 40)
    log_message("STEP 2: Store raw counts in layer")
    log_message("=" * 40)
    adata = store_raw_counts(adata)

    # Step 3: Filter genes
    log_message("\n" + "=" * 40)
    log_message("STEP 3: Filter genes (min_cells=3)")
    log_message("=" * 40)
    adata = filter_genes(adata, min_cells=3)

    # Step 4: Normalize and log transform
    log_message("\n" + "=" * 40)
    log_message("STEP 4: Normalize and log transform")
    log_message("=" * 40)
    adata = normalize_and_log(adata, target_sum=10000)

    # Step 5: Select highly variable genes
    log_message("\n" + "=" * 40)
    log_message("STEP 5: Select highly variable genes")
    log_message("=" * 40)
    adata = select_highly_variable_genes(adata, n_top_genes=3000)

    # Step 6: Scale data
    log_message("\n" + "=" * 40)
    log_message("STEP 6: Scale data")
    log_message("=" * 40)
    adata = scale_data(adata, max_value=10)

    # Step 7: PCA
    log_message("\n" + "=" * 40)
    log_message("STEP 7: PCA")
    log_message("=" * 40)
    adata = run_pca(adata, n_comps=50)

    # Step 8: Cell cycle scoring
    log_message("\n" + "=" * 40)
    log_message("STEP 8: Cell cycle scoring")
    log_message("=" * 40)
    adata = score_cell_cycle(adata, S_GENES, G2M_GENES)

    # Step 9: Compute neighbors
    log_message("\n" + "=" * 40)
    log_message("STEP 9: Compute neighbors")
    log_message("=" * 40)
    adata = compute_neighbors(adata, n_neighbors=15, n_pcs=30)

    # Step 10: UMAP
    log_message("\n" + "=" * 40)
    log_message("STEP 10: UMAP")
    log_message("=" * 40)
    adata = run_umap(adata)

    # Step 11: Leiden clustering
    log_message("\n" + "=" * 40)
    log_message("STEP 11: Leiden clustering")
    log_message("=" * 40)
    adata = run_leiden_clustering(adata, resolution=0.5)

    # Step 12: Save UMAP plots
    log_message("\n" + "=" * 40)
    log_message("STEP 12: Save UMAP plots")
    log_message("=" * 40)
    save_umap_plots(adata)

    # Step 13: Plot epithelial markers
    log_message("\n" + "=" * 40)
    log_message("STEP 13: Plot epithelial markers")
    log_message("=" * 40)
    plot_epithelial_markers(adata, EPITHELIAL_MARKERS)

    # Step 14: Save preprocessed data
    log_message("\n" + "=" * 40)
    log_message("STEP 14: Save preprocessed data")
    log_message("=" * 40)
    save_preprocessed_data(adata)

    # Final summary
    log_message("\n" + "=" * 60)
    log_message("PREPROCESSING PIPELINE COMPLETE")
    log_message("=" * 60)
    log_message(f"Final data: {adata.n_obs} cells, {adata.n_vars} genes")
    log_message(f"HVGs: {adata.var['highly_variable'].sum()}")
    log_message(f"Clusters: {adata.obs['leiden'].nunique()}")
    log_message(f"\nOutputs:")
    log_message(f"  Preprocessed data: {PREPROCESSED_FILE}")
    log_message(f"  UMAP figures: {UMAP_FIGURES_DIR}/")

    return adata


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}", flush=True)
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
