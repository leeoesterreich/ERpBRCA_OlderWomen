#!/usr/bin/env python3
"""
03_generate_oncoplot.py

Generate oncoplot visualization from filtered variant data.

This script:
1. Loads filtered variant data from 01_parse_vep.py output
2. Maps rat gene IDs to human gene symbols via pybiomart
3. Filters for cancer-related genes using brca_genelist.csv
4. Generates a publication-quality oncoplot figure

Author: Alexander Chang
Date: 2024-07
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap
from pybiomart import Dataset

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"
FIGURES_DIR = SCRIPT_DIR / "figures"

# Sample age groups
YOUNG_SAMPLES = {"157", "158", "167"}
OLD_SAMPLES = {"102", "107", "116"}


# ============================================================================
# Helper Functions
# ============================================================================


def chunk_list(lst, chunk_size):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i : i + chunk_size]


def find_homologs(gene_ids):
    """Query rat-to-human homologs from Ensembl via pybiomart."""
    try:
        rat_dataset = Dataset(
            name="rnorvegicus_gene_ensembl", host="http://www.ensembl.org"
        )
        result = rat_dataset.query(
            attributes=["ensembl_gene_id", "hsapiens_homolog_ensembl_gene"],
            filters={"link_ensembl_gene_id": gene_ids},
        )
        if not result.empty:
            return result
    except Exception as e:
        print(f"Error finding homologs: {e}")
    return pd.DataFrame()


def find_symbols(ensembl_ids):
    """Query human gene symbols from Ensembl IDs via pybiomart."""
    try:
        dataset = Dataset(name="hsapiens_gene_ensembl", host="http://www.ensembl.org")
        result = dataset.query(
            attributes=["ensembl_gene_id", "external_gene_name"],
            filters={"link_ensembl_gene_id": ensembl_ids},
        )
        if not result.empty:
            return result
    except Exception as e:
        print(f"Error finding gene symbols: {e}")
    return pd.DataFrame()


def map_to_human_symbols(df):
    """Map rat gene IDs to human gene symbols.

    Returns DataFrame with columns: rat_gene_id, human_gene_id, gene_symbol
    """
    gene_ids = df["Gene"].unique().tolist()
    print(f"Mapping {len(gene_ids)} unique rat gene IDs...")

    # Step 1: Get human homologs
    all_homologs = pd.DataFrame()
    for chunk in chunk_list(gene_ids, 200):
        result = find_homologs(chunk)
        all_homologs = pd.concat([all_homologs, result], ignore_index=True)

    if all_homologs.empty:
        return pd.DataFrame()

    all_homologs = all_homologs.rename(
        columns={
            "Gene stable ID": "rat_gene_id",
            "Human gene stable ID": "human_gene_id",
        }
    )
    all_homologs = all_homologs.dropna(subset=["human_gene_id"])
    print(f"Found {len(all_homologs)} rat-to-human homolog mappings")

    # Step 2: Get human gene symbols
    human_ids = all_homologs["human_gene_id"].unique().tolist()
    all_symbols = pd.DataFrame()
    for chunk in chunk_list(human_ids, 200):
        result = find_symbols(chunk)
        all_symbols = pd.concat([all_symbols, result], ignore_index=True)

    if all_symbols.empty:
        return all_homologs

    all_symbols = all_symbols.rename(
        columns={"Gene stable ID": "human_gene_id", "Gene name": "gene_symbol"}
    )
    all_symbols = all_symbols.dropna(subset=["gene_symbol"])

    # Merge to get final mapping
    mapping = all_homologs.merge(all_symbols, on="human_gene_id", how="inner")
    print(f"Final mapping: {len(mapping)} gene symbol mappings")

    return mapping


def filter_cancer_genes(df, gene_list_file):
    """Filter for genes in cancer gene list."""
    cancer_genes = pd.read_csv(gene_list_file)
    cancer_gene_set = set(cancer_genes["Gene"].tolist())
    return df[df["gene_symbol"].isin(cancer_gene_set)]


def generate_oncoplot_figure(df):
    """Generate publication-quality oncoplot figure."""
    FIGURES_DIR.mkdir(exist_ok=True)

    # Extract consequence and impact
    df = df.copy()
    df["IMPACT"] = df["Extra"].str.extract(r"IMPACT=([^;]+)")

    # Priority: HIGH > MODERATE
    impact_priority = {"HIGH": 1, "MODERATE": 2}
    df["IMPACT_PRIORITY"] = df["IMPACT"].map(impact_priority).fillna(3)
    df = df.sort_values(["Sample", "gene_symbol", "IMPACT_PRIORITY"])

    # Aggregate to one consequence per gene per sample
    agg_df = df.groupby(["Sample", "gene_symbol"]).first().reset_index()

    # Pivot to matrix format
    plot_data = agg_df.pivot_table(
        values="Consequence", index="gene_symbol", columns="Sample", aggfunc="first"
    )

    # Sort genes by mutation frequency
    gene_freq = plot_data.notna().sum(axis=1).sort_values(ascending=False)
    plot_data = plot_data.reindex(gene_freq.index)

    # Sort samples: Old first, then Young
    sample_order = sorted(
        plot_data.columns, key=lambda x: (x[:3] in YOUNG_SAMPLES, x)
    )
    plot_data = plot_data[sample_order]

    # Create consequence-to-number mapping
    unique_consequences = df["Consequence"].dropna().unique()
    consequence_map = {cons: i + 1 for i, cons in enumerate(unique_consequences)}

    # Convert to numeric matrix (0 for missing)
    plot_numeric = plot_data.map(
        lambda x: consequence_map.get(x, 0) if pd.notna(x) else 0
    )

    # Create colormap with white for missing
    n_colors = len(unique_consequences) + 1
    colors = plt.cm.tab20(np.linspace(0, 1, n_colors))
    colors[0] = [1, 1, 1, 1]  # White for missing
    custom_cmap = ListedColormap(colors)

    # Plot
    num_genes = plot_numeric.shape[0]
    fig_height = max(num_genes * 0.4, 6)

    fig, ax = plt.subplots(figsize=(10, fig_height))
    sns.heatmap(plot_numeric, cmap=custom_cmap, cbar=False, linewidths=0.5, ax=ax)

    # Add gridlines
    lines = []
    for i in range(plot_numeric.shape[0] + 1):
        lines.append(((0, i), (plot_numeric.shape[1], i)))
    for j in range(plot_numeric.shape[1] + 1):
        lines.append(((j, 0), (j, plot_numeric.shape[0])))
    line_segments = LineCollection(lines, color="gray", linewidths=0.5, alpha=0.5)
    ax.add_collection(line_segments)

    ax.set_title("Oncoplot", fontsize=16)
    ax.set_xlabel("Samples", fontsize=14)
    ax.set_ylabel("Genes", fontsize=14)

    # Relabel x-axis with age suffix
    new_labels = []
    for label in ax.get_xticklabels():
        sample_id = label.get_text()[:3]
        suffix = "_Y" if sample_id in YOUNG_SAMPLES else "_O"
        new_labels.append(sample_id + suffix)
    ax.set_xticklabels(new_labels, fontsize=12, rotation=45, ha="right")
    ax.tick_params(axis="y", labelsize=10)

    # Legend
    legend_elements = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=colors[consequence_map[cons]],
            edgecolor="none",
            label=cons.replace("_", " ").title(),
        )
        for cons in unique_consequences
    ]
    ax.legend(
        handles=legend_elements,
        title="Consequence",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=9,
    )

    plt.tight_layout()

    fig.savefig(FIGURES_DIR / "oncoplot.svg", format="svg", bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "oncoplot.png", format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved oncoplot figure: {num_genes} genes x {plot_data.shape[1]} samples")

    return plot_data


# ============================================================================
# Main
# ============================================================================


def main():
    print("=== Generating Oncoplot Data and Figure ===")

    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)

    # Load filtered variant data
    combined_file = OUTPUT_DIR / "all_samples_filtered.csv"
    if not combined_file.exists():
        raise FileNotFoundError("Run 01_parse_vep.py first")

    df = pd.read_csv(combined_file)
    print(f"Loaded {len(df)} filtered variants")

    # Map to human gene symbols
    print("Mapping rat genes to human symbols...")
    mapping = map_to_human_symbols(df)

    if mapping.empty:
        raise ValueError("No gene mappings found - check pybiomart connection")

    # Merge mapping with variant data
    df = df.merge(mapping, left_on="Gene", right_on="rat_gene_id", how="inner")
    print(f"Mapped {len(df)} variants to human symbols")

    # Filter for cancer genes
    gene_list_file = DATA_DIR / "brca_genelist.csv"
    if not gene_list_file.exists():
        raise FileNotFoundError(f"Gene list not found: {gene_list_file}")

    df_cancer = filter_cancer_genes(df, gene_list_file)
    print(
        f"Filtered to {len(df_cancer)} variants in {df_cancer['gene_symbol'].nunique()} cancer genes"
    )

    if df_cancer.empty:
        print("Warning: No cancer genes found in data")
        return

    # Generate figure
    plot_data = generate_oncoplot_figure(df_cancer)

    # Save data outputs
    plot_data.to_csv(OUTPUT_DIR / "oncoplot_data.csv")

    summary = pd.DataFrame(
        {
            "gene": plot_data.index,
            "n_samples_mutated": plot_data.notna().sum(axis=1),
            "mutation_frequency": plot_data.notna().sum(axis=1) / plot_data.shape[1],
        }
    )
    summary.to_csv(OUTPUT_DIR / "oncoplot_summary.csv", index=False)

    print(f"Saved oncoplot data to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
