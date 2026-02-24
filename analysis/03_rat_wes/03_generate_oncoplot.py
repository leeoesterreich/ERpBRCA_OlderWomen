#!/usr/bin/env python3
"""Generate oncoplot data from filtered VEP output."""

import pandas as pd
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent / "outputs"


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
