#!/usr/bin/env python
# analysis/07_organoid_single_cell/06_cell_cycle.py
# Cell cycle distribution analysis with chi-square tests
#
# Inputs:
#   - data/processed/final.h5ad (from 05_pathways.py)
#
# Outputs:
#   - outputs/cell_cycle/cell_cycle_proportions.csv
#   - outputs/cell_cycle/cell_cycle_chi_square.csv
#   - figures/cell_cycle/phase_stacked_bar.png/.pdf
#   - figures/cell_cycle/cell_cycle_scores_boxplot.png/.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _figure_config import save_fig
import seaborn as sns

import scanpy as sc
import pandas as pd
import numpy as np
import os
from datetime import datetime
from scipy.stats import chi2_contingency

from _config import (
    OUTPUT_DIR, FIGURES_DIR, h5ad_path, check_file_exists, standardize_treatments
)

# =============================================================================
# Configuration
# =============================================================================
# Input path
FINAL_FILE = h5ad_path("final.h5ad")

# Results and figures directories
CELL_CYCLE_RESULTS_DIR = OUTPUT_DIR / "cell_cycle"
CELL_CYCLE_FIGURES_DIR = FIGURES_DIR / "cell_cycle"

# Treatment order for consistent plotting
TREATMENT_ORDER = ["Vehicle", "E1", "E1+fulv", "E1+HSD17B7i", "E2", "E2+fulv", "E2+HSD17B7i"]

# Phase colors (G1=green, S=yellow, G2M=red)
PHASE_COLORS = ["#2ecc71", "#f1c40f", "#e74c3c"]
PHASE_ORDER = ["G1", "S", "G2M"]

# Key comparisons for chi-square tests
COMPARISONS = [
    ("E1", "E2"),
    ("E1", "E1+HSD17B7i"),
    ("E1+HSD17B7i", "E2+HSD17B7i"),
]

# Scanpy settings
sc.settings.verbosity = 2


def log_message(msg):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {msg}")


def load_final_data():
    """Load final AnnData object.

    Returns
    -------
    adata : AnnData
        Final AnnData object with cell cycle scores
    """
    log_message(f"Loading final data from {FINAL_FILE}")
    check_file_exists(FINAL_FILE, "Final file")

    adata = sc.read_h5ad(FINAL_FILE)
    standardize_treatments(adata)
    log_message(f"  Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Check for required columns
    required_cols = ["treatment", "phase", "S_score", "G2M_score"]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in adata.obs: {missing_cols}")

    log_message(f"  Found required columns: {required_cols}")

    return adata


def calculate_phase_distribution(adata):
    """Calculate phase counts and proportions per treatment.

    Parameters
    ----------
    adata : AnnData
        AnnData with treatment and phase columns

    Returns
    -------
    phase_counts : DataFrame
        Phase counts per treatment (rows=treatments, cols=phases)
    phase_props : DataFrame
        Phase proportions per treatment (rows=treatments, cols=phases)
    """
    log_message("\nCalculating phase distribution per treatment...")

    # Phase counts per treatment
    phase_counts = adata.obs.groupby(["treatment", "phase"]).size().unstack(fill_value=0)

    # Ensure all phases are present
    for phase in PHASE_ORDER:
        if phase not in phase_counts.columns:
            phase_counts[phase] = 0
    phase_counts = phase_counts[PHASE_ORDER]

    # Calculate proportions
    phase_props = phase_counts.div(phase_counts.sum(axis=1), axis=0)

    # Reorder rows by treatment order (only include treatments present in data)
    available_treatments = [t for t in TREATMENT_ORDER if t in phase_counts.index]
    phase_counts = phase_counts.loc[available_treatments]
    phase_props = phase_props.loc[available_treatments]

    log_message("  Phase counts per treatment:")
    for treatment in available_treatments:
        counts_str = ", ".join([f"{p}={phase_counts.loc[treatment, p]}" for p in PHASE_ORDER])
        log_message(f"    {treatment}: {counts_str}")

    log_message("\n  Phase proportions per treatment:")
    for treatment in available_treatments:
        props_str = ", ".join([f"{p}={phase_props.loc[treatment, p]:.1%}" for p in PHASE_ORDER])
        log_message(f"    {treatment}: {props_str}")

    return phase_counts, phase_props


def chi2_test(phase_counts, treatment1, treatment2):
    """Perform chi-square test comparing phase distributions between two treatments.

    Parameters
    ----------
    phase_counts : DataFrame
        Phase counts per treatment
    treatment1 : str
        First treatment name
    treatment2 : str
        Second treatment name

    Returns
    -------
    chi2 : float
        Chi-square statistic
    p : float
        P-value
    dof : int
        Degrees of freedom
    expected : ndarray
        Expected frequencies
    """
    if treatment1 not in phase_counts.index:
        raise ValueError(f"Treatment '{treatment1}' not found in data")
    if treatment2 not in phase_counts.index:
        raise ValueError(f"Treatment '{treatment2}' not found in data")

    t1_counts = phase_counts.loc[treatment1]
    t2_counts = phase_counts.loc[treatment2]
    contingency = pd.DataFrame([t1_counts, t2_counts])

    chi2, p, dof, expected = chi2_contingency(contingency)

    return chi2, p, dof, expected


def run_chi_square_tests(phase_counts):
    """Run chi-square tests for key comparisons.

    Parameters
    ----------
    phase_counts : DataFrame
        Phase counts per treatment

    Returns
    -------
    results : DataFrame
        Chi-square test results for all comparisons
    """
    log_message("\nRunning chi-square tests for key comparisons...")

    results = []

    for treatment1, treatment2 in COMPARISONS:
        # Check if treatments are available
        if treatment1 not in phase_counts.index:
            log_message(f"  Skipping {treatment1} vs {treatment2}: {treatment1} not found")
            continue
        if treatment2 not in phase_counts.index:
            log_message(f"  Skipping {treatment1} vs {treatment2}: {treatment2} not found")
            continue

        chi2, p, dof, expected = chi2_test(phase_counts, treatment1, treatment2)

        # Determine significance
        if p < 0.001:
            sig = "***"
        elif p < 0.01:
            sig = "**"
        elif p < 0.05:
            sig = "*"
        else:
            sig = "ns"

        results.append({
            "Comparison": f"{treatment1} vs {treatment2}",
            "Chi2": chi2,
            "p-value": p,
            "dof": dof,
            "Significance": sig
        })

        log_message(f"  {treatment1} vs {treatment2}: chi2={chi2:.2f}, p={p:.4e} {sig}")

    results_df = pd.DataFrame(results)
    return results_df


def plot_stacked_bar(phase_props, phase_counts):
    """Create stacked bar plot of phase proportions by treatment.

    Parameters
    ----------
    phase_props : DataFrame
        Phase proportions per treatment
    phase_counts : DataFrame
        Phase counts per treatment
    """
    log_message("\nCreating stacked bar plot of phase proportions...")

    fig, ax = plt.subplots(figsize=(12, 6))

    # Create stacked bar plot
    bottom = np.zeros(len(phase_props))
    x = np.arange(len(phase_props))

    for i, (phase, color) in enumerate(zip(PHASE_ORDER, PHASE_COLORS)):
        values = phase_props[phase].values
        bars = ax.bar(
            x,
            values,
            bottom=bottom,
            label=phase,
            color=color,
            edgecolor="white",
            linewidth=0.5
        )
        bottom += values

    # Customize plot
    ax.set_xlabel("Treatment", fontsize=12)
    ax.set_ylabel("Proportion", fontsize=12)
    ax.set_title("Cell Cycle Phase Distribution by Treatment", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(phase_props.index, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.legend(title="Phase", loc="upper right", bbox_to_anchor=(1.15, 1))

    # Add cell count annotations below x-axis
    for i, treatment in enumerate(phase_props.index):
        n_cells = phase_counts.loc[treatment].sum()
        ax.annotate(
            f"n={n_cells}",
            xy=(i, -0.08),
            xycoords=("data", "axes fraction"),
            ha="center",
            va="top",
            fontsize=10,
            color="gray"
        )

    out_path_png = CELL_CYCLE_FIGURES_DIR / "phase_stacked_bar.png"
    out_path_pdf = CELL_CYCLE_FIGURES_DIR / "phase_stacked_bar.pdf"
    save_fig(out_path_png)
    save_fig(out_path_pdf)
    plt.close()
    log_message(f"  Saved: {out_path_png}")


def plot_cell_cycle_scores(adata):
    """Create boxplots for S_score and G2M_score by treatment.

    Parameters
    ----------
    adata : AnnData
        AnnData with S_score and G2M_score in obs
    """
    log_message("\nCreating boxplots for cell cycle scores...")

    # Get available treatments in the correct order
    available_treatments = [t for t in TREATMENT_ORDER if t in adata.obs["treatment"].unique()]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot S_score
    ax = axes[0]
    plot_df = pd.DataFrame({
        "S Score": adata.obs["S_score"].values,
        "Treatment": adata.obs["treatment"].values
    })
    sns.boxplot(
        data=plot_df,
        x="Treatment",
        y="S Score",
        hue="Treatment",
        order=available_treatments,
        ax=ax,
        palette="Set2",
        showfliers=False,
        legend=False
    )
    sns.stripplot(
        data=plot_df,
        x="Treatment",
        y="S Score",
        order=available_treatments,
        ax=ax,
        color="black",
        alpha=0.3,
        size=2
    )
    ax.set_xlabel("")
    ax.set_ylabel("S Score", fontsize=12)
    ax.set_title("S Phase Score by Treatment", fontsize=14)
    ax.tick_params(axis="x", rotation=45)

    # Plot G2M_score
    ax = axes[1]
    plot_df = pd.DataFrame({
        "G2M Score": adata.obs["G2M_score"].values,
        "Treatment": adata.obs["treatment"].values
    })
    sns.boxplot(
        data=plot_df,
        x="Treatment",
        y="G2M Score",
        hue="Treatment",
        order=available_treatments,
        ax=ax,
        palette="Set2",
        showfliers=False,
        legend=False
    )
    sns.stripplot(
        data=plot_df,
        x="Treatment",
        y="G2M Score",
        order=available_treatments,
        ax=ax,
        color="black",
        alpha=0.3,
        size=2
    )
    ax.set_xlabel("")
    ax.set_ylabel("G2M Score", fontsize=12)
    ax.set_title("G2/M Phase Score by Treatment", fontsize=14)
    ax.tick_params(axis="x", rotation=45)

    out_path_png = CELL_CYCLE_FIGURES_DIR / "cell_cycle_scores_boxplot.png"
    out_path_pdf = CELL_CYCLE_FIGURES_DIR / "cell_cycle_scores_boxplot.pdf"
    save_fig(out_path_png)
    save_fig(out_path_pdf)
    plt.close()
    log_message(f"  Saved: {out_path_png}")


def save_proportions(phase_counts, phase_props, chi2_results):
    """Save cell cycle proportions and chi-square results to CSV.

    Parameters
    ----------
    phase_counts : DataFrame
        Phase counts per treatment
    phase_props : DataFrame
        Phase proportions per treatment
    chi2_results : DataFrame
        Chi-square test results
    """
    log_message("\nSaving cell cycle proportions and statistics...")

    # Combine counts and proportions into a single output
    output_df = phase_counts.copy()
    output_df.columns = [f"{col}_count" for col in output_df.columns]

    props_renamed = phase_props.copy()
    props_renamed.columns = [f"{col}_proportion" for col in props_renamed.columns]

    output_df = pd.concat([output_df, props_renamed], axis=1)
    output_df["total_cells"] = phase_counts.sum(axis=1)
    output_df.index.name = "treatment"

    # Save proportions
    props_path = CELL_CYCLE_RESULTS_DIR / "cell_cycle_proportions.csv"
    output_df.to_csv(props_path)
    log_message(f"  Saved: {props_path}")

    # Save chi-square results
    if len(chi2_results) > 0:
        chi2_path = CELL_CYCLE_RESULTS_DIR / "cell_cycle_chi_square.csv"
        chi2_results.to_csv(chi2_path, index=False)
        log_message(f"  Saved: {chi2_path}")


def main():
    """Main cell cycle analysis pipeline."""
    log_message("=" * 60)
    log_message("Cell Cycle Distribution Analysis")
    log_message("=" * 60)

    # Ensure output directories exist
    CELL_CYCLE_FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    CELL_CYCLE_RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Load final data
    log_message("\n" + "=" * 40)
    log_message("STEP 1: Load final data")
    log_message("=" * 40)
    adata = load_final_data()

    # Step 2: Calculate phase distribution
    log_message("\n" + "=" * 40)
    log_message("STEP 2: Calculate phase distribution")
    log_message("=" * 40)
    phase_counts, phase_props = calculate_phase_distribution(adata)

    # Step 3: Run chi-square tests
    log_message("\n" + "=" * 40)
    log_message("STEP 3: Run chi-square tests")
    log_message("=" * 40)
    chi2_results = run_chi_square_tests(phase_counts)

    # Step 4: Create stacked bar plot
    log_message("\n" + "=" * 40)
    log_message("STEP 4: Create stacked bar plot")
    log_message("=" * 40)
    plot_stacked_bar(phase_props, phase_counts)

    # Step 5: Create cell cycle score boxplots
    log_message("\n" + "=" * 40)
    log_message("STEP 5: Create cell cycle score boxplots")
    log_message("=" * 40)
    plot_cell_cycle_scores(adata)

    # Step 6: Save results
    log_message("\n" + "=" * 40)
    log_message("STEP 6: Save results")
    log_message("=" * 40)
    save_proportions(phase_counts, phase_props, chi2_results)

    # Final summary
    log_message("\n" + "=" * 60)
    log_message("CELL CYCLE ANALYSIS COMPLETE")
    log_message("=" * 60)

    log_message("\nSummary:")
    log_message(f"  Total cells: {adata.n_obs}")
    log_message(f"  Treatments analyzed: {len(phase_counts)}")
    log_message(f"  Chi-square comparisons: {len(chi2_results)}")

    # Log significant comparisons
    if len(chi2_results) > 0:
        sig_results = chi2_results[chi2_results["Significance"] != "ns"]
        if len(sig_results) > 0:
            log_message("\n  Significant comparisons (p<0.05):")
            for _, row in sig_results.iterrows():
                log_message(f"    {row['Comparison']}: p={row['p-value']:.4e} {row['Significance']}")
        else:
            log_message("\n  No significant differences in phase distribution")

    log_message(f"\nOutputs:")
    log_message(f"  Figures: {CELL_CYCLE_FIGURES_DIR}/")
    log_message(f"  Proportions: {CELL_CYCLE_RESULTS_DIR}/cell_cycle_proportions.csv")
    log_message(f"  Chi-square tests: {CELL_CYCLE_RESULTS_DIR}/cell_cycle_chi_square.csv")

    return adata, phase_counts, phase_props, chi2_results


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}", flush=True)
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
