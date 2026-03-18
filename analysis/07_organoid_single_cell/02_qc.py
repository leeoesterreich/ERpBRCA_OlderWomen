#!/usr/bin/env python
# analysis/07_organoid_single_cell/02_qc.py
# Quality control for organoid scRNA-seq data (Cell Ranger Flex output)
#
# Inputs:
#   - Cell Ranger per_sample_outs for OS01-OS07 (via SAMPLE_POOL_MAP)
#   - configs/sample_metadata.csv
#
# Outputs:
#   - outputs/raw_merged.h5ad
#   - outputs/qc_filtered.h5ad
#   - outputs/qc/cell_counts.csv
#   - figures/qc/qc_violin_pre_filter.png/.pdf
#   - figures/qc/qc_violin_post_filter.png/.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _figure_config import save_fig

import scanpy as sc
import scrublet as scr
import pandas as pd
import numpy as np
import os
from datetime import datetime

from _config import (
    OUTPUT_DIR, FIGURES_DIR, CONFIGS_DIR, DATA_PATHS,
    SAMPLE_IDS, SAMPLE_POOL_MAP, MIN_GENES, MIN_COUNTS, MAX_PCT_MITO,
    MIN_CELLS_FOR_SCRUBLET, h5ad_path, check_file_exists
)

# =============================================================================
# Configuration
# =============================================================================
METADATA_FILE = CONFIGS_DIR / "sample_metadata.csv"

# Output directories
QC_RESULTS_DIR = OUTPUT_DIR / "qc"
QC_FIGURES_DIR = FIGURES_DIR / "qc"

# Output files
RAW_MERGED_FILE = h5ad_path("raw_merged.h5ad")
QC_FILTERED_FILE = h5ad_path("qc_filtered.h5ad")
CELL_COUNTS_FILE = QC_RESULTS_DIR / "cell_counts.csv"

# Scanpy settings
sc.settings.verbosity = 2
sc.settings.figdir = str(QC_FIGURES_DIR)
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)


def log_message(msg):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}")


def load_sample(sample_id):
    """Load Cell Ranger output for a single sample.

    Parameters
    ----------
    sample_id : str
        Sample identifier (e.g., 'OS01')

    Returns
    -------
    adata : AnnData
        AnnData object with sample_id added to obs
    """
    pool_out = SAMPLE_POOL_MAP[sample_id]
    matrix_path = f"{pool_out}/{sample_id}/count/sample_filtered_feature_bc_matrix.h5"

    if not os.path.exists(matrix_path):
        raise FileNotFoundError(f"Matrix file not found: {matrix_path}")

    log_message(f"Loading {sample_id} from {matrix_path}")
    adata = sc.read_10x_h5(matrix_path)

    # Make var_names unique (in case of duplicates)
    adata.var_names_make_unique()

    # Add sample_id to observations
    adata.obs['sample_id'] = sample_id

    # Make cell barcodes unique by prepending sample_id
    adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]

    log_message(f"  {sample_id}: {adata.n_obs} cells, {adata.n_vars} genes")

    return adata


def load_and_merge_samples():
    """Load all samples and merge into single AnnData.

    Returns
    -------
    adata : AnnData
        Merged AnnData object
    """
    log_message("Loading and merging samples...")

    adatas = []
    for sample_id in SAMPLE_IDS:
        try:
            adata = load_sample(sample_id)
            adatas.append(adata)
        except FileNotFoundError as e:
            log_message(f"WARNING: {e}")
            continue

    if len(adatas) == 0:
        raise RuntimeError("No samples could be loaded!")

    log_message(f"Merging {len(adatas)} samples...")
    adata = sc.concat(adatas, join='outer')

    log_message(f"Merged data: {adata.n_obs} cells, {adata.n_vars} genes")

    return adata


def add_sample_metadata(adata, metadata_file):
    """Add sample metadata to AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData object with sample_id in obs
    metadata_file : str
        Path to sample metadata CSV

    Returns
    -------
    adata : AnnData
        AnnData with metadata added to obs
    """
    log_message(f"Adding sample metadata from {metadata_file}")

    metadata = pd.read_csv(metadata_file)
    metadata = metadata.set_index('sample_id')

    # Map metadata columns to cells based on sample_id
    for col in metadata.columns:
        adata.obs[col] = adata.obs['sample_id'].map(metadata[col])

    log_message(f"  Added columns: {list(metadata.columns)}")

    return adata


def calculate_qc_metrics(adata):
    """Calculate QC metrics including mitochondrial percentage.

    Parameters
    ----------
    adata : AnnData
        AnnData object

    Returns
    -------
    adata : AnnData
        AnnData with QC metrics in obs
    """
    log_message("Calculating QC metrics...")

    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    n_mt_genes = adata.var['mt'].sum()
    log_message(f"  Found {n_mt_genes} mitochondrial genes")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    # Log summary statistics
    log_message("  QC metric summaries:")
    log_message(f"    n_genes_by_counts: median={adata.obs['n_genes_by_counts'].median():.0f}, "
                f"range=[{adata.obs['n_genes_by_counts'].min():.0f}, {adata.obs['n_genes_by_counts'].max():.0f}]")
    log_message(f"    total_counts: median={adata.obs['total_counts'].median():.0f}, "
                f"range=[{adata.obs['total_counts'].min():.0f}, {adata.obs['total_counts'].max():.0f}]")
    log_message(f"    pct_counts_mt: median={adata.obs['pct_counts_mt'].median():.1f}%, "
                f"range=[{adata.obs['pct_counts_mt'].min():.1f}%, {adata.obs['pct_counts_mt'].max():.1f}%]")

    return adata


def plot_qc_violin(adata, suffix, title_prefix=""):
    """Create QC violin plots.

    Parameters
    ----------
    adata : AnnData
        AnnData object with QC metrics
    suffix : str
        Suffix for output filename (e.g., 'pre_filter', 'post_filter')
    title_prefix : str
        Prefix for plot titles
    """
    log_message(f"Creating QC violin plots ({suffix})...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # n_genes_by_counts
    sc.pl.violin(
        adata,
        keys='n_genes_by_counts',
        groupby='sample_id',
        ax=axes[0],
        show=False,
        rotation=45
    )
    axes[0].set_title(f'{title_prefix}Genes per Cell')
    axes[0].set_ylabel('n_genes_by_counts')
    if suffix == 'pre_filter':
        axes[0].axhline(y=MIN_GENES, color='red', linestyle='--', alpha=0.7, label=f'threshold={MIN_GENES}')
        axes[0].legend()

    # total_counts
    sc.pl.violin(
        adata,
        keys='total_counts',
        groupby='sample_id',
        ax=axes[1],
        show=False,
        rotation=45
    )
    axes[1].set_title(f'{title_prefix}Total Counts per Cell')
    axes[1].set_ylabel('total_counts')
    if suffix == 'pre_filter':
        axes[1].axhline(y=MIN_COUNTS, color='red', linestyle='--', alpha=0.7, label=f'threshold={MIN_COUNTS}')
        axes[1].legend()

    # pct_counts_mt
    sc.pl.violin(
        adata,
        keys='pct_counts_mt',
        groupby='sample_id',
        ax=axes[2],
        show=False,
        rotation=45
    )
    axes[2].set_title(f'{title_prefix}Mitochondrial %')
    axes[2].set_ylabel('pct_counts_mt')
    if suffix == 'pre_filter':
        axes[2].axhline(y=MAX_PCT_MITO, color='red', linestyle='--', alpha=0.7, label=f'threshold={MAX_PCT_MITO}%')
        axes[2].legend()

    # Save figure (PNG + PDF)
    out_path_png = QC_FIGURES_DIR / f"qc_violin_{suffix}.png"
    out_path_pdf = QC_FIGURES_DIR / f"qc_violin_{suffix}.pdf"
    save_fig(out_path_png)
    save_fig(out_path_pdf)
    plt.close()
    log_message(f"  Saved: {out_path_png}")
    log_message(f"  Saved: {out_path_pdf}")


def apply_qc_filters(adata):
    """Apply QC filters to remove low-quality cells.

    Parameters
    ----------
    adata : AnnData
        AnnData object with QC metrics

    Returns
    -------
    adata : AnnData
        Filtered AnnData object
    """
    log_message("Applying QC filters...")
    n_cells_before = adata.n_obs

    # Filter cells with too few genes
    adata = adata[adata.obs['n_genes_by_counts'] >= MIN_GENES, :].copy()
    n_after_genes = adata.n_obs
    log_message(f"  After min_genes >= {MIN_GENES}: {n_after_genes} cells "
                f"(removed {n_cells_before - n_after_genes})")

    # Filter cells with too few counts
    adata = adata[adata.obs['total_counts'] >= MIN_COUNTS, :].copy()
    n_after_counts = adata.n_obs
    log_message(f"  After min_counts >= {MIN_COUNTS}: {n_after_counts} cells "
                f"(removed {n_after_genes - n_after_counts})")

    # Filter cells with high mitochondrial percentage
    adata = adata[adata.obs['pct_counts_mt'] < MAX_PCT_MITO, :].copy()
    n_after_mito = adata.n_obs
    log_message(f"  After pct_mito < {MAX_PCT_MITO}%: {n_after_mito} cells "
                f"(removed {n_after_counts - n_after_mito})")

    log_message(f"  Total removed by QC filters: {n_cells_before - adata.n_obs} cells "
                f"({100 * (n_cells_before - adata.n_obs) / n_cells_before:.1f}%)")

    return adata


def run_scrublet_per_sample(adata):
    """Run scrublet doublet detection for each sample separately.

    Parameters
    ----------
    adata : AnnData
        AnnData object with sample_id in obs

    Returns
    -------
    adata : AnnData
        AnnData with doublet_score and predicted_doublet in obs
    """
    log_message("Running scrublet doublet detection per sample...")

    # Initialize columns
    adata.obs['doublet_score'] = 0.0
    adata.obs['predicted_doublet'] = False

    for sample_id in adata.obs['sample_id'].unique():
        sample_mask = adata.obs['sample_id'] == sample_id
        n_cells = sample_mask.sum()

        log_message(f"  Processing {sample_id} ({n_cells} cells)...")

        # Skip samples with too few cells
        if n_cells < MIN_CELLS_FOR_SCRUBLET:
            log_message(f"    Skipping scrublet (< {MIN_CELLS_FOR_SCRUBLET} cells)")
            continue

        # Get count matrix for this sample
        sample_adata = adata[sample_mask, :].copy()
        counts_matrix = sample_adata.X

        try:
            # Initialize scrublet
            scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

            # Run scrublet
            doublet_scores, predicted_doublets = scrub.scrub_doublets(
                min_counts=2,
                min_cells=3,
                min_gene_variability_pctl=85,
                n_prin_comps=30,
                verbose=False
            )

            # Store results
            adata.obs.loc[sample_mask, 'doublet_score'] = doublet_scores
            adata.obs.loc[sample_mask, 'predicted_doublet'] = predicted_doublets

            n_doublets = predicted_doublets.sum()
            log_message(f"    Found {n_doublets} doublets ({100 * n_doublets / n_cells:.1f}%)")

        except Exception as e:
            log_message(f"    WARNING: Scrublet failed for {sample_id}: {e}")
            log_message(f"    Setting all cells as singlets for this sample")
            continue

    return adata


def remove_doublets(adata):
    """Remove predicted doublets from the dataset.

    Parameters
    ----------
    adata : AnnData
        AnnData with predicted_doublet in obs

    Returns
    -------
    adata : AnnData
        AnnData with doublets removed
    """
    log_message("Removing predicted doublets...")
    n_before = adata.n_obs
    n_doublets = adata.obs['predicted_doublet'].sum()

    adata = adata[~adata.obs['predicted_doublet'], :].copy()

    log_message(f"  Removed {n_doublets} doublets ({100 * n_doublets / n_before:.1f}%)")
    log_message(f"  Remaining cells: {adata.n_obs}")

    return adata


def create_cell_count_summary(adata_raw, adata_filtered):
    """Create summary of cell counts before and after filtering.

    Parameters
    ----------
    adata_raw : AnnData
        Raw merged AnnData (before QC filtering)
    adata_filtered : AnnData
        Filtered AnnData (after QC and doublet removal)

    Returns
    -------
    summary_df : DataFrame
        Summary of cell counts per sample
    """
    log_message("Creating cell count summary...")

    # Count cells per sample
    raw_counts = adata_raw.obs['sample_id'].value_counts().sort_index()
    filtered_counts = adata_filtered.obs['sample_id'].value_counts().sort_index()

    # Create summary DataFrame
    summary_df = pd.DataFrame({
        'sample_id': raw_counts.index,
        'raw_cells': raw_counts.values,
        'filtered_cells': [filtered_counts.get(s, 0) for s in raw_counts.index]
    })

    summary_df['removed_cells'] = summary_df['raw_cells'] - summary_df['filtered_cells']
    summary_df['pct_retained'] = 100 * summary_df['filtered_cells'] / summary_df['raw_cells']

    # Add totals row
    totals = pd.DataFrame({
        'sample_id': ['TOTAL'],
        'raw_cells': [summary_df['raw_cells'].sum()],
        'filtered_cells': [summary_df['filtered_cells'].sum()],
        'removed_cells': [summary_df['removed_cells'].sum()],
        'pct_retained': [100 * summary_df['filtered_cells'].sum() / summary_df['raw_cells'].sum()]
    })
    summary_df = pd.concat([summary_df, totals], ignore_index=True)

    return summary_df


def main():
    """Main QC pipeline."""
    log_message("=" * 60)
    log_message("Starting QC pipeline for organoid scRNA-seq")
    log_message("=" * 60)

    # Ensure output directories exist
    QC_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    QC_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Load and merge samples
    log_message("\n" + "=" * 40)
    log_message("STEP 1: Load and merge samples")
    log_message("=" * 40)
    adata = load_and_merge_samples()

    # Step 2: Add sample metadata
    log_message("\n" + "=" * 40)
    log_message("STEP 2: Add sample metadata")
    log_message("=" * 40)
    adata = add_sample_metadata(adata, str(METADATA_FILE))

    # Step 3: Calculate QC metrics
    log_message("\n" + "=" * 40)
    log_message("STEP 3: Calculate QC metrics")
    log_message("=" * 40)
    adata = calculate_qc_metrics(adata)

    # Step 4: Save raw merged data
    log_message("\n" + "=" * 40)
    log_message("STEP 4: Save raw merged data")
    log_message("=" * 40)
    log_message(f"Saving raw merged data to {RAW_MERGED_FILE}")
    adata.write_h5ad(RAW_MERGED_FILE)
    log_message(f"  Saved: {adata.n_obs} cells, {adata.n_vars} genes")

    # Keep a copy for cell count summary
    adata_raw = adata.copy()

    # Step 5: Pre-filter QC violin plots
    log_message("\n" + "=" * 40)
    log_message("STEP 5: Pre-filter QC plots")
    log_message("=" * 40)
    plot_qc_violin(adata, 'pre_filter', title_prefix='Pre-filter: ')

    # Step 6: Apply QC filters
    log_message("\n" + "=" * 40)
    log_message("STEP 6: Apply QC filters")
    log_message("=" * 40)
    adata = apply_qc_filters(adata)

    # Step 7: Run scrublet per sample
    log_message("\n" + "=" * 40)
    log_message("STEP 7: Scrublet doublet detection")
    log_message("=" * 40)
    adata = run_scrublet_per_sample(adata)

    # Step 8: Remove doublets
    log_message("\n" + "=" * 40)
    log_message("STEP 8: Remove doublets")
    log_message("=" * 40)
    adata = remove_doublets(adata)

    # Step 9: Post-filter QC violin plots
    log_message("\n" + "=" * 40)
    log_message("STEP 9: Post-filter QC plots")
    log_message("=" * 40)
    plot_qc_violin(adata, 'post_filter', title_prefix='Post-filter: ')

    # Step 10: Save cell count summary
    log_message("\n" + "=" * 40)
    log_message("STEP 10: Save cell count summary")
    log_message("=" * 40)
    summary_df = create_cell_count_summary(adata_raw, adata)
    summary_df.to_csv(CELL_COUNTS_FILE, index=False)
    log_message(f"Saved cell count summary to {CELL_COUNTS_FILE}")
    log_message("\nCell count summary:")
    print(summary_df.to_string(index=False))

    # Step 11: Save QC filtered data
    log_message("\n" + "=" * 40)
    log_message("STEP 11: Save QC filtered data")
    log_message("=" * 40)
    log_message(f"Saving QC filtered data to {QC_FILTERED_FILE}")
    adata.write_h5ad(QC_FILTERED_FILE)
    log_message(f"  Saved: {adata.n_obs} cells, {adata.n_vars} genes")

    # Final summary
    log_message("\n" + "=" * 60)
    log_message("QC PIPELINE COMPLETE")
    log_message("=" * 60)
    log_message(f"Raw cells: {adata_raw.n_obs}")
    log_message(f"Filtered cells: {adata.n_obs}")
    log_message(f"Retention rate: {100 * adata.n_obs / adata_raw.n_obs:.1f}%")
    log_message(f"\nOutputs:")
    log_message(f"  Raw merged data: {RAW_MERGED_FILE}")
    log_message(f"  QC filtered data: {QC_FILTERED_FILE}")
    log_message(f"  Cell counts: {CELL_COUNTS_FILE}")
    log_message(f"  QC figures: {QC_FIGURES_DIR}/")

    return adata


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}", flush=True)
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
