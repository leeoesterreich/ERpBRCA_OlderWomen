#!/usr/bin/env python
# analysis/07_organoid_single_cell/05_pathways.py
# GSEA and decoupler pathway/TF activity analysis
#
# Inputs:
#   - data/processed/preprocessed.h5ad (from 03_preprocess.py)
#   - outputs/de/*_deseq2_results.csv (from 04_pseudobulk_de.py)
#
# Outputs:
#   - outputs/pathways/gsea_<comparison>.csv
#   - figures/pathways/tf_boxplots.png/.pdf
#   - figures/pathways/progeny_heatmap.png/.pdf
#   - figures/pathways/gsea_summary_heatmap.png/.pdf
#   - data/processed/final.h5ad

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import pandas as pd
import numpy as np
import os
from datetime import datetime
from glob import glob

import gseapy as gp
import decoupler as dc

from _config import (
    OUTPUT_DIR, FIGURES_DIR, h5ad_path, check_file_exists
)

# =============================================================================
# Configuration
# =============================================================================
# Input/Output paths
PREPROCESSED_FILE = h5ad_path("preprocessed.h5ad")
FINAL_FILE = h5ad_path("final.h5ad")

# DE results directory
DE_RESULTS_DIR = OUTPUT_DIR / "de"

# Results and figures directories
PATHWAYS_RESULTS_DIR = OUTPUT_DIR / "pathways"
PATHWAYS_FIGURES_DIR = FIGURES_DIR / "pathways"

# GSEA parameters
GSEA_GENE_SETS = "MSigDB_Hallmark_2020"
GSEA_MIN_SIZE = 15
GSEA_MAX_SIZE = 500
GSEA_PERMUTATION_NUM = 1000
GSEA_SEED = 42

# TFs to plot
TFS_TO_PLOT = ["ESR1", "E2F1", "E2F4", "MYC", "TP53"]

# Comparisons (same as in DE script)
COMPARISONS = [
    "E1_vs_E2",
    "E1_vs_E1_HSD17B7i",
    "E1_HSD17B7i_vs_E2_HSD17B7i",
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
        Preprocessed AnnData object
    """
    log_message(f"Loading preprocessed data from {PREPROCESSED_FILE}")
    check_file_exists(PREPROCESSED_FILE, "Preprocessed file")

    adata = sc.read_h5ad(PREPROCESSED_FILE)
    log_message(f"  Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    return adata


def load_de_results():
    """Load DE results from CSV files.

    Returns
    -------
    de_results : dict
        Dictionary mapping comparison name to DE results DataFrame
    """
    log_message("Loading DE results...")

    de_results = {}
    de_files = sorted(DE_RESULTS_DIR.glob("*_deseq2_results.csv"))

    if not de_files:
        log_message(f"  WARNING: No DE result files found in {DE_RESULTS_DIR}")
        return de_results

    for de_file in de_files:
        # Extract comparison name from filename
        comparison_name = de_file.stem.replace("_deseq2_results", "")

        log_message(f"  Loading: {de_file.name}")
        de_df = pd.read_csv(de_file, index_col=0)
        de_results[comparison_name] = de_df
        log_message(f"    {len(de_df)} genes")

    return de_results


def run_preranked_gsea(de_results):
    """Run pre-ranked GSEA for all DE comparisons.

    Parameters
    ----------
    de_results : dict
        Dictionary mapping comparison name to DE results DataFrame

    Returns
    -------
    gsea_results : dict
        Dictionary mapping comparison name to GSEA results DataFrame
    """
    log_message("\nRunning pre-ranked GSEA...")
    log_message(f"  Gene sets: {GSEA_GENE_SETS}")
    log_message(f"  Min size: {GSEA_MIN_SIZE}, Max size: {GSEA_MAX_SIZE}")
    log_message(f"  Permutations: {GSEA_PERMUTATION_NUM}")

    gsea_results = {}

    for comparison_name, de_df in de_results.items():
        log_message(f"\n  Processing: {comparison_name}")

        # Create ranking: log2FoldChange * -log10(pvalue)
        # Handle missing or zero p-values
        de_df = de_df.dropna(subset=["log2FoldChange", "pvalue"])
        de_df = de_df[de_df["pvalue"] > 0]  # Remove zero p-values

        if len(de_df) == 0:
            log_message(f"    WARNING: No valid genes for GSEA in {comparison_name}")
            continue

        # Calculate ranking metric
        de_df["rank"] = de_df["log2FoldChange"] * -np.log10(de_df["pvalue"].clip(lower=1e-300))

        # Sort by rank
        ranking = de_df["rank"].sort_values(ascending=False)
        log_message(f"    Using {len(ranking)} genes for ranking")

        try:
            # Run pre-ranked GSEA
            gsea_res = gp.prerank(
                rnk=ranking,
                gene_sets=GSEA_GENE_SETS,
                outdir=None,
                min_size=GSEA_MIN_SIZE,
                max_size=GSEA_MAX_SIZE,
                permutation_num=GSEA_PERMUTATION_NUM,
                seed=GSEA_SEED,
                verbose=False
            )

            results_df = gsea_res.res2d
            gsea_results[comparison_name] = results_df

            # Log significant pathways
            sig_pathways = results_df[results_df["FDR q-val"] < 0.25]
            n_up = len(sig_pathways[sig_pathways["NES"] > 0])
            n_down = len(sig_pathways[sig_pathways["NES"] < 0])
            log_message(f"    Significant pathways (FDR<0.25): {len(sig_pathways)} ({n_up} up, {n_down} down)")

            # Save results
            out_path = PATHWAYS_RESULTS_DIR / f"gsea_{comparison_name}.csv"
            results_df.to_csv(out_path)
            log_message(f"    Saved: {out_path}")

        except Exception as e:
            log_message(f"    ERROR running GSEA: {e}")
            import traceback
            traceback.print_exc()

    return gsea_results


def run_decoupler_progeny(adata):
    """Run PROGENy pathway activity analysis using decoupler.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData object

    Returns
    -------
    adata : AnnData
        AnnData with PROGENy scores in obsm['progeny']
    """
    log_message("\nRunning PROGENy pathway activity analysis...")

    # Get PROGENy network (decoupler 2.x API)
    log_message("  Fetching PROGENy network (human, top 300)...")
    progeny = dc.op.progeny(organism="human", top=300)
    log_message(f"    Network: {len(progeny)} interactions")
    log_message(f"    Pathways: {progeny['source'].nunique()}")

    # Run MLM (Multivariate Linear Model)
    log_message("  Running MLM...")
    dc.mt.mlm(data=adata, net=progeny, raw=False)

    # Store results
    for key in ["score_mlm", "mlm_estimate", "mlm_scores"]:
        if key in adata.obsm:
            adata.obsm["progeny"] = adata.obsm[key].copy()
            break
    log_message(f"  Stored PROGENy scores in adata.obsm['progeny']")
    log_message(f"    Shape: {adata.obsm['progeny'].shape}")

    return adata


def run_decoupler_dorothea(adata):
    """Run DoRothEA TF activity analysis using decoupler.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData object

    Returns
    -------
    adata : AnnData
        AnnData with DoRothEA scores in obsm['dorothea']
    """
    log_message("\nRunning DoRothEA TF activity analysis...")

    # Get DoRothEA network (decoupler 2.x API)
    log_message("  Fetching DoRothEA network (human, levels A,B,C)...")
    dorothea = dc.op.dorothea(organism="human", levels=["A", "B", "C"])
    log_message(f"    Network: {len(dorothea)} interactions")
    log_message(f"    TFs: {dorothea['source'].nunique()}")

    # Run MLM (Multivariate Linear Model)
    log_message("  Running MLM...")
    dc.mt.mlm(data=adata, net=dorothea, raw=False)

    # Store results
    for key in ["score_mlm", "mlm_estimate", "mlm_scores"]:
        if key in adata.obsm:
            adata.obsm["dorothea"] = adata.obsm[key].copy()
            break
    log_message(f"  Stored DoRothEA scores in adata.obsm['dorothea']")
    log_message(f"    Shape: {adata.obsm['dorothea'].shape}")

    return adata


def plot_tf_boxplots(adata, tfs_to_plot):
    """Plot TF activity boxplots by treatment.

    Parameters
    ----------
    adata : AnnData
        AnnData with DoRothEA scores in obsm['dorothea']
    tfs_to_plot : list
        List of TF names to plot
    """
    log_message("\nPlotting TF activity boxplots...")

    if "dorothea" not in adata.obsm:
        log_message("  WARNING: DoRothEA scores not found in adata.obsm")
        return

    # Get DoRothEA scores
    tf_scores = adata.obsm["dorothea"]

    # Check which TFs are available
    available_tfs = [tf for tf in tfs_to_plot if tf in tf_scores.columns]
    missing_tfs = [tf for tf in tfs_to_plot if tf not in tf_scores.columns]

    if missing_tfs:
        log_message(f"  WARNING: Missing TFs: {', '.join(missing_tfs)}")

    if not available_tfs:
        log_message("  ERROR: No TFs found, skipping plot")
        return

    log_message(f"  Plotting TFs: {', '.join(available_tfs)}")

    # Create figure
    n_tfs = len(available_tfs)
    ncols = min(3, n_tfs)
    nrows = (n_tfs + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))

    # Handle single TF case
    if n_tfs == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    # Get treatment order for consistent plotting
    treatment_order = sorted(adata.obs["treatment"].unique())

    for i, tf in enumerate(available_tfs):
        ax = axes[i]

        # Create DataFrame for plotting
        plot_df = pd.DataFrame({
            "TF Activity": tf_scores[tf].values,
            "Treatment": adata.obs["treatment"].values
        })

        # Create boxplot
        sns.boxplot(
            data=plot_df,
            x="Treatment",
            y="TF Activity",
            hue="Treatment",
            order=treatment_order,
            ax=ax,
            palette="Set2",
            legend=False
        )

        ax.set_title(f"{tf} Activity", fontsize=12)
        ax.set_xlabel("")
        ax.tick_params(axis='x', rotation=45)

    # Remove empty axes
    for j in range(len(available_tfs), len(axes)):
        axes[j].axis('off')

    plt.tight_layout()
    out_path_png = PATHWAYS_FIGURES_DIR / "tf_boxplots.png"
    out_path_pdf = PATHWAYS_FIGURES_DIR / "tf_boxplots.pdf"
    plt.savefig(out_path_png, dpi=150, bbox_inches='tight')
    plt.savefig(out_path_pdf, bbox_inches='tight')
    plt.close()
    log_message(f"  Saved: {out_path_png}")


def plot_progeny_heatmap(adata):
    """Plot PROGENy pathway activity heatmap (mean per treatment).

    Parameters
    ----------
    adata : AnnData
        AnnData with PROGENy scores in obsm['progeny']
    """
    log_message("\nPlotting PROGENy pathway heatmap...")

    if "progeny" not in adata.obsm:
        log_message("  WARNING: PROGENy scores not found in adata.obsm")
        return

    # Get PROGENy scores
    pathway_scores = adata.obsm["progeny"]

    # Calculate mean per treatment
    treatments = adata.obs["treatment"].unique()
    mean_scores = pd.DataFrame(index=treatments, columns=pathway_scores.columns)

    for treatment in treatments:
        mask = adata.obs["treatment"] == treatment
        mean_scores.loc[treatment] = pathway_scores[mask].mean(axis=0)

    mean_scores = mean_scores.astype(float)

    log_message(f"  Pathways: {mean_scores.shape[1]}")
    log_message(f"  Treatments: {', '.join(treatments)}")

    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 6))

    sns.heatmap(
        mean_scores.T,  # Pathways as rows, treatments as columns
        cmap="RdBu_r",
        center=0,
        annot=True,
        fmt=".2f",
        ax=ax,
        cbar_kws={"label": "Mean Activity Score"}
    )

    ax.set_title("PROGENy Pathway Activities by Treatment", fontsize=14)
    ax.set_xlabel("Treatment", fontsize=12)
    ax.set_ylabel("Pathway", fontsize=12)

    plt.tight_layout()
    out_path_png = PATHWAYS_FIGURES_DIR / "progeny_heatmap.png"
    out_path_pdf = PATHWAYS_FIGURES_DIR / "progeny_heatmap.pdf"
    plt.savefig(out_path_png, dpi=150, bbox_inches='tight')
    plt.savefig(out_path_pdf, bbox_inches='tight')
    plt.close()
    log_message(f"  Saved: {out_path_png}")


def plot_gsea_summary(gsea_results):
    """Plot summary of GSEA results across comparisons.

    Parameters
    ----------
    gsea_results : dict
        Dictionary mapping comparison name to GSEA results DataFrame
    """
    log_message("\nPlotting GSEA summary...")

    if not gsea_results:
        log_message("  WARNING: No GSEA results to plot")
        return

    # Combine results from all comparisons
    all_results = []
    for comparison_name, results_df in gsea_results.items():
        df = results_df.copy()
        df["Comparison"] = comparison_name
        all_results.append(df)

    combined = pd.concat(all_results, ignore_index=True)

    # Get top pathways by significance across all comparisons
    # Use absolute NES for ranking
    combined["abs_NES"] = combined["NES"].abs()
    top_pathways = combined.groupby("Term")["abs_NES"].max().nlargest(15).index.tolist()

    # Filter to top pathways
    plot_data = combined[combined["Term"].isin(top_pathways)]

    # Pivot for heatmap
    pivot_nes = plot_data.pivot(index="Term", columns="Comparison", values="NES")

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 12))

    # Create heatmap with significance markers
    sns.heatmap(
        pivot_nes,
        cmap="RdBu_r",
        center=0,
        annot=True,
        fmt=".2f",
        ax=ax,
        cbar_kws={"label": "Normalized Enrichment Score (NES)"}
    )

    ax.set_title("Top GSEA Pathways Across Comparisons", fontsize=14)
    ax.set_xlabel("Comparison", fontsize=12)
    ax.set_ylabel("Pathway", fontsize=12)

    plt.tight_layout()
    out_path_png = PATHWAYS_FIGURES_DIR / "gsea_summary_heatmap.png"
    out_path_pdf = PATHWAYS_FIGURES_DIR / "gsea_summary_heatmap.pdf"
    plt.savefig(out_path_png, dpi=150, bbox_inches='tight')
    plt.savefig(out_path_pdf, bbox_inches='tight')
    plt.close()
    log_message(f"  Saved: {out_path_png}")


def save_final_data(adata):
    """Save final AnnData with pathway scores.

    Parameters
    ----------
    adata : AnnData
        AnnData with pathway scores
    """
    log_message(f"\nSaving final data to {FINAL_FILE}")
    adata.write_h5ad(FINAL_FILE)
    log_message(f"  Saved: {adata.n_obs} cells, {adata.n_vars} genes")
    log_message(f"  obsm keys: {list(adata.obsm.keys())}")


def main():
    """Main pathway analysis pipeline."""
    log_message("=" * 60)
    log_message("Pathway Analysis Pipeline")
    log_message("=" * 60)

    # Ensure output directories exist
    PATHWAYS_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PATHWAYS_FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Load preprocessed data
    log_message("\n" + "=" * 40)
    log_message("STEP 1: Load preprocessed data")
    log_message("=" * 40)
    adata = load_preprocessed_data()

    # Step 2: Load DE results
    log_message("\n" + "=" * 40)
    log_message("STEP 2: Load DE results")
    log_message("=" * 40)
    de_results = load_de_results()

    # Step 3: Run pre-ranked GSEA
    log_message("\n" + "=" * 40)
    log_message("STEP 3: Run pre-ranked GSEA")
    log_message("=" * 40)
    gsea_results = {}
    if de_results:
        gsea_results = run_preranked_gsea(de_results)
    else:
        log_message("  Skipping GSEA - no DE results available")

    # Step 4: Run PROGENy pathway analysis
    log_message("\n" + "=" * 40)
    log_message("STEP 4: Run PROGENy pathway analysis")
    log_message("=" * 40)
    adata = run_decoupler_progeny(adata)

    # Step 5: Run DoRothEA TF analysis
    log_message("\n" + "=" * 40)
    log_message("STEP 5: Run DoRothEA TF analysis")
    log_message("=" * 40)
    adata = run_decoupler_dorothea(adata)

    # Step 6: Plot TF boxplots
    log_message("\n" + "=" * 40)
    log_message("STEP 6: Plot TF activity boxplots")
    log_message("=" * 40)
    plot_tf_boxplots(adata, TFS_TO_PLOT)

    # Step 7: Plot PROGENy heatmap
    log_message("\n" + "=" * 40)
    log_message("STEP 7: Plot PROGENy heatmap")
    log_message("=" * 40)
    plot_progeny_heatmap(adata)

    # Step 8: Plot GSEA summary
    log_message("\n" + "=" * 40)
    log_message("STEP 8: Plot GSEA summary")
    log_message("=" * 40)
    if gsea_results:
        plot_gsea_summary(gsea_results)
    else:
        log_message("  Skipping GSEA summary plot - no results available")

    # Step 9: Save final data
    log_message("\n" + "=" * 40)
    log_message("STEP 9: Save final data")
    log_message("=" * 40)
    save_final_data(adata)

    # Final summary
    log_message("\n" + "=" * 60)
    log_message("PATHWAY ANALYSIS COMPLETE")
    log_message("=" * 60)

    log_message("\nSummary:")
    log_message(f"  DE comparisons processed: {len(de_results)}")
    log_message(f"  GSEA results generated: {len(gsea_results)}")
    log_message(f"  PROGENy pathways: {adata.obsm['progeny'].shape[1] if 'progeny' in adata.obsm else 0}")
    log_message(f"  DoRothEA TFs: {adata.obsm['dorothea'].shape[1] if 'dorothea' in adata.obsm else 0}")

    log_message(f"\nOutputs:")
    log_message(f"  GSEA results: {PATHWAYS_RESULTS_DIR}/gsea_*.csv")
    log_message(f"  Pathway figures: {PATHWAYS_FIGURES_DIR}/")
    log_message(f"  Final data: {FINAL_FILE}")

    return adata, gsea_results


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}", flush=True)
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
