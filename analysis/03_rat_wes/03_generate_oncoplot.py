#!/usr/bin/env python3
"""Generate oncoplot data and figure from filtered VEP output."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.collections import LineCollection
from pathlib import Path
from pybiomart import Dataset

OUTPUT_DIR = Path(__file__).parent / "outputs"
FIGURES_DIR = Path(__file__).parent / "figures"
DATA_DIR = Path(__file__).parent / "data"

# Sample age groups (confirmed: lower IDs are older rats)
YOUNG_SAMPLES = ['157', '158', '167']
OLD_SAMPLES = ['102', '107', '116']


def extract_sample_id(sample_name: str) -> str:
    """Extract numeric sample ID from various naming formats.

    Handles: '157_tumor_vs_157_spleen', '102', 'Rat_O_102', etc.
    """
    import re
    match = re.search(r'(102|107|116|157|158|167)', sample_name)
    if match:
        return match.group(1)
    return sample_name[:3]  # Fallback to first 3 chars


def chunk_list(lst, chunk_size=200):
    """Yield successive chunks from list."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]


def find_homologs(gene_ids):
    """Map rat Ensembl gene IDs to human homologs."""
    try:
        rat_dataset = Dataset(name='rnorvegicus_gene_ensembl', host='http://www.ensembl.org')
        result = rat_dataset.query(
            attributes=['ensembl_gene_id', 'hsapiens_homolog_ensembl_gene'],
            filters={'link_ensembl_gene_id': gene_ids}
        )
        if not result.empty:
            return result
    except Exception as e:
        print(f"Error finding homologs: {e}")
    return pd.DataFrame()


def find_symbols(ensembl_ids):
    """Map human Ensembl IDs to gene symbols."""
    try:
        dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
        result = dataset.query(
            attributes=['ensembl_gene_id', 'external_gene_name'],
            filters={'link_ensembl_gene_id': ensembl_ids}
        )
        if not result.empty:
            return result
    except Exception as e:
        print(f"Error finding gene symbols: {e}")
    return pd.DataFrame()


def map_to_human_symbols(df):
    """Map rat genes to human gene symbols via homologs."""
    gene_ids = df['Gene'].unique().tolist()

    # Check for cached results
    cache_file = OUTPUT_DIR / "homolog_cache.csv"
    if cache_file.exists():
        print("Loading cached homolog mapping...")
        cache = pd.read_csv(cache_file)
        cached_genes = set(cache['rat_gene_id'].tolist())
        uncached = [g for g in gene_ids if g not in cached_genes]
        if not uncached:
            return cache
        gene_ids = uncached
        print(f"Querying {len(gene_ids)} uncached genes...")
    else:
        cache = pd.DataFrame()

    # Query homologs in chunks
    all_homologs = pd.DataFrame()
    for chunk in chunk_list(gene_ids, 200):
        result = find_homologs(chunk)
        if not result.empty:
            all_homologs = pd.concat([all_homologs, result], ignore_index=True)

    if all_homologs.empty:
        print("Warning: No homologs found for uncached genes")
        # Return cached data if available
        if not cache.empty:
            return cache
        return pd.DataFrame()

    # Rename columns
    all_homologs.columns = ['rat_gene_id', 'human_ensembl_id']
    all_homologs = all_homologs.dropna()

    # Query gene symbols in chunks
    human_ids = all_homologs['human_ensembl_id'].unique().tolist()
    all_symbols = pd.DataFrame()
    for chunk in chunk_list(human_ids, 200):
        result = find_symbols(chunk)
        if not result.empty:
            all_symbols = pd.concat([all_symbols, result], ignore_index=True)

    if all_symbols.empty:
        print("Warning: No gene symbols found for uncached genes")
        # Return cached data if available
        if not cache.empty:
            return cache
        return all_homologs

    all_symbols.columns = ['human_ensembl_id', 'gene_symbol']

    # Merge to get final mapping
    mapping = all_homologs.merge(all_symbols, on='human_ensembl_id', how='left')
    mapping = mapping.dropna(subset=['gene_symbol'])

    # Update cache
    if not cache.empty:
        mapping = pd.concat([cache, mapping], ignore_index=True).drop_duplicates()
    mapping.to_csv(cache_file, index=False)
    print(f"Cached {len(mapping)} gene mappings")

    return mapping


def filter_cancer_genes(df, gene_list_file):
    """Filter for genes in cancer gene list."""
    cancer_genes = pd.read_csv(gene_list_file)
    cancer_gene_set = set(cancer_genes['Gene'].tolist())
    return df[df['gene_symbol'].isin(cancer_gene_set)]


def generate_oncoplot_figure(df):
    """Generate publication-quality oncoplot figure."""
    FIGURES_DIR.mkdir(exist_ok=True)

    # Extract consequence and impact
    df = df.copy()
    # IMPACT is now a separate column (not in Extra)
    if 'IMPACT' not in df.columns:
        df['IMPACT'] = df['Extra'].str.extract(r'IMPACT=([^;]+)')

    # Priority: HIGH > MODERATE
    impact_priority = {'HIGH': 1, 'MODERATE': 2}
    df['IMPACT_PRIORITY'] = df['IMPACT'].map(impact_priority).fillna(3)
    df = df.sort_values(['Sample', 'gene_symbol', 'IMPACT_PRIORITY'])

    # Aggregate to one consequence per gene per sample
    agg_df = df.groupby(['Sample', 'gene_symbol']).first().reset_index()

    # Pivot to matrix format
    plot_data = agg_df.pivot_table(
        values='Consequence',
        index='gene_symbol',
        columns='Sample',
        aggfunc='first'
    )

    # Sort genes by mutation frequency
    gene_freq = plot_data.notna().sum(axis=1).sort_values(ascending=False)
    plot_data = plot_data.reindex(gene_freq.index)

    # Sort samples: Young first (167, 158, 157), then Old (116, 107, 102)
    # Within each group, use DESCENDING numeric order to match manuscript
    sample_order = sorted(plot_data.columns, key=lambda x: (
        extract_sample_id(x) not in YOUNG_SAMPLES,  # True (Old) sorts after False (Young)
        -int(extract_sample_id(x))  # DESCENDING numeric ID within group
    ))
    plot_data = plot_data[sample_order]

    # Create consequence-to-number mapping
    unique_consequences = df['Consequence'].dropna().unique()
    consequence_map = {cons: i+1 for i, cons in enumerate(unique_consequences)}

    # Convert to numeric matrix (0 for missing)
    # Use map for pandas 2.x compatibility
    plot_numeric = plot_data.applymap(lambda x: consequence_map.get(x, 0) if pd.notna(x) else 0)

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
    line_segments = LineCollection(lines, color='gray', linewidths=0.5, alpha=0.5)
    ax.add_collection(line_segments)

    ax.set_title('Oncoplot', fontsize=16)
    ax.set_xlabel('Samples', fontsize=14)
    ax.set_ylabel('Genes', fontsize=14)

    # Relabel x-axis with age suffix using ID-based lookup (not positional)
    # Per legend: Young (157, 158, 167) = '_Y', Old (102, 107, 116) = '_O'
    new_labels = []
    for label in ax.get_xticklabels():
        sample_id = extract_sample_id(label.get_text())
        suffix = '_Y' if sample_id in YOUNG_SAMPLES else '_O'
        new_labels.append(f"{sample_id}{suffix}")
    ax.set_xticklabels(new_labels, fontsize=12, rotation=45, ha='right')
    ax.tick_params(axis='y', labelsize=10)

    # Legend
    legend_elements = [
        plt.Rectangle((0,0), 1, 1, facecolor=colors[consequence_map[cons]],
                      edgecolor='none', label=cons.replace('_', ' ').title())
        for cons in unique_consequences
    ]
    ax.legend(handles=legend_elements, title='Consequence',
             bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)

    plt.tight_layout()

    fig.savefig(FIGURES_DIR / "oncoplot.svg", format='svg', bbox_inches='tight')
    fig.savefig(FIGURES_DIR / "oncoplot.png", format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved oncoplot figure: {num_genes} genes x {plot_data.shape[1]} samples")

    return plot_data


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
    df = df.merge(mapping, left_on='Gene', right_on='rat_gene_id', how='inner')
    print(f"Mapped {len(df)} variants to human symbols")

    # Filter for cancer genes
    gene_list_file = DATA_DIR / "brca_genelist.csv"
    if not gene_list_file.exists():
        raise FileNotFoundError(f"Gene list not found: {gene_list_file}")

    df_cancer = filter_cancer_genes(df, gene_list_file)
    print(f"Filtered to {len(df_cancer)} variants in {df_cancer['gene_symbol'].nunique()} cancer genes")

    if df_cancer.empty:
        print("Warning: No cancer genes found in data")
        return

    # Generate figure
    plot_data = generate_oncoplot_figure(df_cancer)

    # Save data outputs
    plot_data.to_csv(OUTPUT_DIR / "oncoplot_data.csv")

    summary = pd.DataFrame({
        'gene': plot_data.index,
        'n_samples_mutated': plot_data.notna().sum(axis=1),
        'mutation_frequency': plot_data.notna().sum(axis=1) / plot_data.shape[1]
    })
    summary.to_csv(OUTPUT_DIR / "oncoplot_summary.csv", index=False)

    print(f"Saved oncoplot data to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
