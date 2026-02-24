#!/usr/bin/env python3
"""Parse VEP output files and filter for high/moderate impact variants."""

import os
import pandas as pd
from pathlib import Path

# Input directory
VEP_DIR = Path("/bgfs/alee/LO_LAB/General/Lab_Data/20240628_WES_Rat_Neil/results/variant_calling/mutect2/")
OUTPUT_DIR = Path(__file__).parent / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)


def read_vep_output(file_path):
    """Read VEP output file with custom header parsing."""
    header = None
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') and not line.startswith('##'):
                header = line[1:].strip().split('\t')
                break
    if header is None:
        raise ValueError(f"No VEP header found in {file_path}")
    df = pd.read_csv(file_path, comment='#', sep='\t', names=header, header=None)
    return df


def filter_vep_output(df):
    """Filter for Ensembl genes with HIGH or MODERATE impact."""
    ens_filter = df['Gene'].str.startswith('ENS', na=False)
    high_impact = df['Extra'].str.contains('IMPACT=HIGH', na=False)
    moderate_impact = df['Extra'].str.contains('IMPACT=MODERATE', na=False)
    return df[ens_filter & (high_impact | moderate_impact)]


def main():
    print("=== Parsing VEP Output Files ===")

    if not VEP_DIR.exists():
        raise FileNotFoundError(f"VEP input directory not found: {VEP_DIR}")

    dataframes = {}
    for subfolder in os.listdir(VEP_DIR):
        if subfolder.endswith('spleen'):
            vep_path = VEP_DIR / subfolder / f"{subfolder}.mutect2.txt"
            if vep_path.exists():
                df = read_vep_output(vep_path)
                filtered = filter_vep_output(df)
                dataframes[subfolder] = filtered
                print(f"  {subfolder}: {len(df)} -> {len(filtered)} variants")

    if not dataframes:
        raise ValueError(f"No VEP files ending in 'spleen' found in {VEP_DIR}")

    # Save filtered results
    for name, df in dataframes.items():
        df.to_csv(OUTPUT_DIR / f"{name}_filtered.csv", index=False)

    # Save combined for downstream analysis
    combined = pd.concat(dataframes.values(), keys=dataframes.keys(), names=['Sample'])
    combined = combined.reset_index(level='Sample')
    combined.to_csv(OUTPUT_DIR / "all_samples_filtered.csv", index=False)

    print(f"\nSaved {len(dataframes)} filtered files to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
