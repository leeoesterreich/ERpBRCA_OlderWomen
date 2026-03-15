#!/usr/bin/env python
# analysis/07_organoid_single_cell/04_pseudobulk_de.py
# Pseudobulk differential expression analysis using PyDESeq2
#
# Inputs:
#   - data/processed/preprocessed.h5ad (from 03_preprocess.py)
#
# Outputs:
#   - outputs/de/<comparison>_deseq2_results.csv
#   - figures/de/volcano_<comparison>.png/.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scanpy as sc
import pandas as pd
import numpy as np
import os
from datetime import datetime
from scipy import sparse

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from _config import (
    OUTPUT_DIR, FIGURES_DIR, h5ad_path, check_file_exists
)

# =============================================================================
# Configuration
# =============================================================================
# Input/Output paths
PREPROCESSED_FILE = h5ad_path("preprocessed.h5ad")

# Results and figures directories
DE_RESULTS_DIR = OUTPUT_DIR / "de"
DE_FIGURES_DIR = FIGURES_DIR / "de"

# Differential expression thresholds
PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 0.5
MIN_GENE_COUNT = 10  # Minimum total counts across samples to keep gene

# Comparisons to run: (name, condition1, condition2)
# PyDESeq2 contrast format: log2FC = log2(condition1 / condition2)
# So positive log2FC means higher in condition1
COMPARISONS = [
    ("E1_vs_E2", "E1", "E2"),
    ("E1_vs_E1_HSD17B7i", "E1", "E1+HSD17B7i"),
    ("E1_HSD17B7i_vs_E2_HSD17B7i", "E1+HSD17B7i", "E2+HSD17B7i"),
]

# Scanpy settings
sc.settings.verbosity = 2


def log_message(msg):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}")


def load_preprocessed_data():
    """Load preprocessed AnnData object.

    Returns
    -------
    adata : AnnData
        Preprocessed AnnData with raw counts in layers["counts"]
    """
    log_message(f"Loading preprocessed data from {PREPROCESSED_FILE}")
    check_file_exists(PREPROCESSED_FILE, "Preprocessed file")

    adata = sc.read_h5ad(PREPROCESSED_FILE)
    log_message(f"  Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Verify counts layer exists
    if "counts" not in adata.layers:
        raise ValueError("Raw counts not found in adata.layers['counts']")

    log_message("  Found raw counts in layers['counts']")

    return adata


def create_pseudobulk(adata):
    """Create pseudobulk count matrix by summing counts per sample.

    Parameters
    ----------
    adata : AnnData
        AnnData with raw counts in layers["counts"] and sample_id in obs

    Returns
    -------
    pseudobulk_df : DataFrame
        Pseudobulk counts matrix (samples x genes)
    sample_meta : DataFrame
        Sample-level metadata
    """
    log_message("Creating pseudobulk count matrix...")

    # Get raw counts
    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        counts = counts.toarray()

    # Get sample IDs
    sample_ids = adata.obs["sample_id"].unique()
    log_message(f"  Found {len(sample_ids)} samples: {', '.join(sorted(sample_ids))}")

    # Aggregate counts per sample
    pseudobulk_data = {}
    sample_cell_counts = {}

    for sample_id in sample_ids:
        mask = adata.obs["sample_id"] == sample_id
        sample_counts = counts[mask, :].sum(axis=0)

        # Ensure 1D array
        if hasattr(sample_counts, 'A1'):
            sample_counts = sample_counts.A1
        elif len(sample_counts.shape) > 1:
            sample_counts = sample_counts.flatten()

        pseudobulk_data[sample_id] = sample_counts
        sample_cell_counts[sample_id] = mask.sum()

    # Create DataFrame (samples x genes)
    pseudobulk_df = pd.DataFrame(
        pseudobulk_data,
        index=adata.var_names
    ).T  # Transpose so samples are rows, genes are columns

    log_message(f"  Pseudobulk matrix shape: {pseudobulk_df.shape[0]} samples x {pseudobulk_df.shape[1]} genes")

    # Create sample-level metadata
    sample_meta_list = []
    for sample_id in pseudobulk_df.index:
        # Get metadata from first cell of this sample
        sample_mask = adata.obs["sample_id"] == sample_id
        first_cell_idx = adata.obs.index[sample_mask][0]
        sample_info = adata.obs.loc[first_cell_idx].to_dict()
        sample_info["n_cells"] = sample_cell_counts[sample_id]
        sample_meta_list.append(sample_info)

    sample_meta = pd.DataFrame(sample_meta_list, index=pseudobulk_df.index)

    # Log sample information
    log_message("  Sample summary:")
    for sample_id in sorted(pseudobulk_df.index):
        n_cells = sample_meta.loc[sample_id, "n_cells"]
        treatment = sample_meta.loc[sample_id, "treatment"]
        total_counts = pseudobulk_df.loc[sample_id].sum()
        log_message(f"    {sample_id}: {treatment}, {n_cells} cells, {total_counts:.0f} total counts")

    return pseudobulk_df, sample_meta


def filter_low_expressed_genes(pseudobulk_df, min_count=10):
    """Filter genes with low total counts across samples.

    Parameters
    ----------
    pseudobulk_df : DataFrame
        Pseudobulk counts (samples x genes)
    min_count : int
        Minimum total counts across samples

    Returns
    -------
    pseudobulk_filtered : DataFrame
        Filtered pseudobulk counts
    """
    log_message(f"Filtering genes with total count < {min_count}...")

    gene_totals = pseudobulk_df.sum(axis=0)
    keep_genes = gene_totals >= min_count

    n_before = pseudobulk_df.shape[1]
    pseudobulk_filtered = pseudobulk_df.loc[:, keep_genes]
    n_after = pseudobulk_filtered.shape[1]

    log_message(f"  Kept {n_after}/{n_before} genes ({100 * n_after / n_before:.1f}%)")

    return pseudobulk_filtered


def run_deseq2_comparison(pseudobulk_df, sample_meta, condition1, condition2, comparison_name):
    """Run DESeq2 for a specific comparison.

    Parameters
    ----------
    pseudobulk_df : DataFrame
        Pseudobulk counts (samples x genes)
    sample_meta : DataFrame
        Sample metadata with 'treatment' column
    condition1 : str
        First condition (numerator in fold change)
    condition2 : str
        Second condition (denominator in fold change)
    comparison_name : str
        Name for this comparison

    Returns
    -------
    results_df : DataFrame
        DESeq2 results with log2FoldChange, pvalue, padj, etc.
    """
    log_message(f"\nRunning DESeq2: {comparison_name}")
    log_message(f"  Contrast: {condition1} vs {condition2}")

    # Filter to only samples in this comparison
    samples_to_use = sample_meta["treatment"].isin([condition1, condition2])
    if samples_to_use.sum() == 0:
        log_message(f"  ERROR: No samples found for conditions {condition1} and {condition2}")
        return None

    pseudobulk_subset = pseudobulk_df.loc[samples_to_use].copy()
    meta_subset = sample_meta.loc[samples_to_use].copy()

    # Log samples being used
    log_message(f"  Using {len(pseudobulk_subset)} samples:")
    for sample_id in pseudobulk_subset.index:
        treatment = meta_subset.loc[sample_id, "treatment"]
        log_message(f"    {sample_id}: {treatment}")

    # Check we have at least 2 samples per condition
    condition_counts = meta_subset["treatment"].value_counts()
    for cond in [condition1, condition2]:
        if cond not in condition_counts:
            log_message(f"  ERROR: No samples found for condition '{cond}'")
            return None
        if condition_counts[cond] < 1:
            log_message(f"  WARNING: Only {condition_counts[cond]} sample(s) for {cond}")

    # Ensure counts are integers
    pseudobulk_subset = pseudobulk_subset.astype(int)

    # Filter genes that are all zeros in this subset
    gene_sums = pseudobulk_subset.sum(axis=0)
    nonzero_genes = gene_sums > 0
    pseudobulk_subset = pseudobulk_subset.loc[:, nonzero_genes]
    log_message(f"  Using {pseudobulk_subset.shape[1]} genes with non-zero counts")

    try:
        # Create DESeq2 dataset
        dds = DeseqDataSet(
            counts=pseudobulk_subset,
            metadata=meta_subset[["treatment"]],
            design_factors="treatment"
        )

        # Run DESeq2 pipeline
        log_message("  Running DESeq2 dispersion estimation...")
        dds.deseq2()

        # Get results for the specific contrast
        # PyDESeq2 contrast: ["factor", "condition1", "condition2"]
        # log2FC = log2(condition1 / condition2), so positive = higher in condition1
        log_message("  Extracting results...")
        stat_res = DeseqStats(dds, contrast=["treatment", condition1, condition2])
        stat_res.summary()

        results_df = stat_res.results_df.copy()

        # Add gene names as column (index is gene name)
        results_df["gene"] = results_df.index

        # Sort by adjusted p-value
        results_df = results_df.sort_values("padj")

        # Count significant DEGs
        sig_mask = (results_df["padj"] < PADJ_THRESHOLD) & (results_df["log2FoldChange"].abs() > LOG2FC_THRESHOLD)
        n_sig = sig_mask.sum()
        n_up = ((results_df["padj"] < PADJ_THRESHOLD) &
                (results_df["log2FoldChange"] > LOG2FC_THRESHOLD)).sum()
        n_down = ((results_df["padj"] < PADJ_THRESHOLD) &
                  (results_df["log2FoldChange"] < -LOG2FC_THRESHOLD)).sum()

        log_message(f"  Results: {len(results_df)} genes tested")
        log_message(f"  Significant DEGs (padj<{PADJ_THRESHOLD}, |log2FC|>{LOG2FC_THRESHOLD}): {n_sig}")
        log_message(f"    Up in {condition1}: {n_up}")
        log_message(f"    Down in {condition1} (up in {condition2}): {n_down}")

        return results_df

    except Exception as e:
        log_message(f"  ERROR running DESeq2: {e}")
        import traceback
        traceback.print_exc()
        return None


def volcano_plot(results_df, title, filename_stem, padj_threshold=0.05, log2fc_threshold=0.5):
    """Create volcano plot for DE results.

    Parameters
    ----------
    results_df : DataFrame
        DESeq2 results with log2FoldChange and padj columns
    title : str
        Plot title
    filename_stem : str
        Output filename stem (without extension)
    padj_threshold : float
        Adjusted p-value threshold for significance
    log2fc_threshold : float
        Log2 fold change threshold for significance
    """
    log_message(f"  Creating volcano plot: {filename_stem}")

    # Remove rows with NaN values
    plot_df = results_df.dropna(subset=["log2FoldChange", "padj"]).copy()

    if len(plot_df) == 0:
        log_message("    No valid data for volcano plot")
        return

    # Calculate -log10(padj)
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].clip(lower=1e-300))

    # Classify points
    sig_up = (plot_df["padj"] < padj_threshold) & (plot_df["log2FoldChange"] > log2fc_threshold)
    sig_down = (plot_df["padj"] < padj_threshold) & (plot_df["log2FoldChange"] < -log2fc_threshold)
    non_sig = ~(sig_up | sig_down)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot non-significant points
    ax.scatter(
        plot_df.loc[non_sig, "log2FoldChange"],
        plot_df.loc[non_sig, "neg_log10_padj"],
        c="gray",
        alpha=0.5,
        s=20,
        label=f"Not significant (n={non_sig.sum()})"
    )

    # Plot significant up-regulated
    ax.scatter(
        plot_df.loc[sig_up, "log2FoldChange"],
        plot_df.loc[sig_up, "neg_log10_padj"],
        c="red",
        alpha=0.7,
        s=30,
        label=f"Up (n={sig_up.sum()})"
    )

    # Plot significant down-regulated
    ax.scatter(
        plot_df.loc[sig_down, "log2FoldChange"],
        plot_df.loc[sig_down, "neg_log10_padj"],
        c="blue",
        alpha=0.7,
        s=30,
        label=f"Down (n={sig_down.sum()})"
    )

    # Add threshold lines
    ax.axhline(y=-np.log10(padj_threshold), color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.axvline(x=log2fc_threshold, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.axvline(x=-log2fc_threshold, color="black", linestyle="--", linewidth=1, alpha=0.5)

    # Labels for top genes
    top_genes = plot_df.loc[sig_up | sig_down].nsmallest(10, "padj")
    for idx, row in top_genes.iterrows():
        ax.annotate(
            idx,
            xy=(row["log2FoldChange"], row["neg_log10_padj"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
            alpha=0.8
        )

    # Formatting
    ax.set_xlabel("log2 Fold Change", fontsize=12)
    ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(loc="upper right")

    # Set reasonable axis limits
    xlim = max(abs(plot_df["log2FoldChange"].min()), abs(plot_df["log2FoldChange"].max()))
    xlim = min(xlim, 10)  # Cap at +/- 10
    ax.set_xlim(-xlim - 0.5, xlim + 0.5)

    plt.tight_layout()
    out_path_png = DE_FIGURES_DIR / f"{filename_stem}.png"
    out_path_pdf = DE_FIGURES_DIR / f"{filename_stem}.pdf"
    plt.savefig(out_path_png, dpi=150, bbox_inches="tight")
    plt.savefig(out_path_pdf, bbox_inches="tight")
    plt.close()

    log_message(f"    Saved: {out_path_png}")


def save_results(results_df, comparison_name):
    """Save DE results to CSV file.

    Parameters
    ----------
    results_df : DataFrame
        DESeq2 results
    comparison_name : str
        Name of comparison for filename
    """
    out_path = DE_RESULTS_DIR / f"{comparison_name}_deseq2_results.csv"
    results_df.to_csv(out_path)
    log_message(f"  Saved results: {out_path}")


def main():
    """Main pseudobulk DE pipeline."""
    log_message("=" * 60)
    log_message("Pseudobulk Differential Expression Analysis")
    log_message("=" * 60)

    # Ensure output directories exist
    DE_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    DE_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Load preprocessed data
    log_message("\n" + "=" * 40)
    log_message("STEP 1: Load preprocessed data")
    log_message("=" * 40)
    adata = load_preprocessed_data()

    # Step 2: Create pseudobulk counts
    log_message("\n" + "=" * 40)
    log_message("STEP 2: Create pseudobulk counts")
    log_message("=" * 40)
    pseudobulk_df, sample_meta = create_pseudobulk(adata)

    # Step 3: Filter low-expressed genes
    log_message("\n" + "=" * 40)
    log_message("STEP 3: Filter low-expressed genes")
    log_message("=" * 40)
    pseudobulk_filtered = filter_low_expressed_genes(pseudobulk_df, min_count=MIN_GENE_COUNT)

    # Step 4: Run DESeq2 for each comparison
    log_message("\n" + "=" * 40)
    log_message("STEP 4: Run DESeq2 comparisons")
    log_message("=" * 40)

    all_results = {}
    for comparison_name, condition1, condition2 in COMPARISONS:
        results = run_deseq2_comparison(
            pseudobulk_filtered,
            sample_meta,
            condition1,
            condition2,
            comparison_name
        )

        if results is not None:
            all_results[comparison_name] = results

            # Save results to CSV
            save_results(results, comparison_name)

            # Create volcano plot
            volcano_plot(
                results,
                title=f"Differential Expression: {condition1} vs {condition2}",
                filename_stem=f"volcano_{comparison_name}",
                padj_threshold=PADJ_THRESHOLD,
                log2fc_threshold=LOG2FC_THRESHOLD
            )

    # Final summary
    log_message("\n" + "=" * 60)
    log_message("PSEUDOBULK DE ANALYSIS COMPLETE")
    log_message("=" * 60)

    log_message("\nSummary of significant DEGs (padj<0.05, |log2FC|>0.5):")
    for comparison_name, results in all_results.items():
        sig_mask = (results["padj"] < PADJ_THRESHOLD) & (results["log2FoldChange"].abs() > LOG2FC_THRESHOLD)
        n_sig = sig_mask.sum()
        n_up = ((results["padj"] < PADJ_THRESHOLD) & (results["log2FoldChange"] > LOG2FC_THRESHOLD)).sum()
        n_down = ((results["padj"] < PADJ_THRESHOLD) & (results["log2FoldChange"] < -LOG2FC_THRESHOLD)).sum()
        log_message(f"  {comparison_name}: {n_sig} DEGs ({n_up} up, {n_down} down)")

    log_message(f"\nOutputs:")
    log_message(f"  DE results: {DE_RESULTS_DIR}/")
    log_message(f"  Volcano plots: {DE_FIGURES_DIR}/")

    return all_results


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}", flush=True)
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
