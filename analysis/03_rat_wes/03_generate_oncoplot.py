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

# Sample age groups
YOUNG_SAMPLES = ['102', '107', '116']
OLD_SAMPLES = ['157', '158', '167']


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
        print("Warning: No homologs found")
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
        print("Warning: No gene symbols found")
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


def load_filtered_data():
    """Load all filtered variant data."""
    combined_file = OUTPUT_DIR / "all_samples_filtered.csv"
    if not combined_file.exists():
        raise FileNotFoundError("Run 01_parse_vep.py first")
    return pd.read_csv(combined_file)


def generate_oncoplot_data(df):
    """Generate oncoplot matrix format."""
    # Extract consequence impact
    df['IMPACT'] = df['Extra'].str.extract(r'IMPACT=([^;]+)')
    df['Consequence'] = df['Consequence'].fillna('Unknown')

    # Create pivot table: genes x samples
    impact_priority = {'HIGH': 1, 'MODERATE': 2}
    df['IMPACT_PRIORITY'] = df['IMPACT'].map(impact_priority).fillna(3)

    df_sorted = df.sort_values(['Sample', 'Gene', 'IMPACT_PRIORITY'])
    agg_df = df_sorted.groupby(['Sample', 'Gene']).first().reset_index()

    # Pivot to matrix format
    oncoplot_matrix = agg_df.pivot_table(
        values='Consequence',
        index='Gene',
        columns='Sample',
        aggfunc='first'
    )

    # Calculate mutation frequency and sort
    freq = oncoplot_matrix.notna().sum(axis=1).sort_values(ascending=False)
    oncoplot_matrix = oncoplot_matrix.reindex(freq.index)

    return oncoplot_matrix


def main():
    print("=== Generating Oncoplot Data ===")

    df = load_filtered_data()
    print(f"Loaded {len(df)} filtered variants")

    oncoplot = generate_oncoplot_data(df)

    # Save to CSV
    oncoplot.to_csv(OUTPUT_DIR / "oncoplot_data.csv")
    print(f"Saved oncoplot data: {oncoplot.shape[0]} genes x {oncoplot.shape[1]} samples")

    # Also save summary statistics
    summary = pd.DataFrame({
        'gene': oncoplot.index,
        'n_samples_mutated': oncoplot.notna().sum(axis=1),
        'mutation_frequency': oncoplot.notna().sum(axis=1) / oncoplot.shape[1]
    })
    summary.to_csv(OUTPUT_DIR / "oncoplot_summary.csv", index=False)


if __name__ == "__main__":
    main()
