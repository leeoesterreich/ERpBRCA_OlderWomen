#!/usr/bin/env python3
"""Run SigProfiler to extract COSMIC mutational signatures."""

import os
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Sample age groups (confirmed: lower IDs are older rats)
YOUNG_SAMPLES = ['157', '158', '167']
OLD_SAMPLES = ['102', '107', '116']

OUTPUT_DIR = Path(__file__).parent / "outputs"
FIGURES_DIR = Path(__file__).parent / "figures"

# SigProfiler paths
VCF_INPUT = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/RatWES_Mutect2_VCF_Input")
SIGPROFILER_OUTPUT = OUTPUT_DIR / "sigprofiler"


def run_sigprofiler():
    """Run SigProfiler COSMIC signature assignment."""
    from SigProfilerAssignment import Analyzer as Analyze

    print("=== Running SigProfiler Assignment ===")

    Analyze.cosmic_fit(
        samples=str(VCF_INPUT),
        output=str(SIGPROFILER_OUTPUT),
        input_type="vcf",
        context_type="96",
        genome_build="rn6",
        cosmic_version=3.4
    )


def extract_sample_id(sample_name: str) -> str:
    """Extract numeric sample ID from various naming formats.

    Handles: 'Rat_O_102', 'Rat_Y_157', '102_tumor_vs_102_spleen', '102', etc.
    """
    import re
    # Look for 3-digit number (102, 107, 116, 157, 158, 167)
    match = re.search(r'(102|107|116|157|158|167)', sample_name)
    if match:
        return match.group(1)
    return sample_name


def get_age_group(sample_id: str) -> str:
    """Get biological age group for a sample ID.

    Confirmed by user:
    - 102, 107, 116 = Old rats
    - 157, 158, 167 = Young rats
    """
    if sample_id in ['102', '107', '116']:
        return 'Old'
    elif sample_id in ['157', '158', '167']:
        return 'Young'
    return 'Unknown'


def parse_signatures():
    """Parse SigProfiler output to CSV format."""
    activities_file = SIGPROFILER_OUTPUT / "Assignment_Solution/Activities/Assignment_Solution_Activities.txt"

    if not activities_file.exists():
        raise FileNotFoundError(f"SigProfiler output not found: {activities_file}")

    df = pd.read_csv(activities_file, sep='\t')
    df = df.set_index('Samples')

    # Drop columns where all values are zero
    df = df.loc[:, (df != 0).any(axis=0)]

    # Add age group annotation using ID-based lookup (not positional)
    df['sample_id'] = df.index.map(extract_sample_id)
    df['age_group'] = df['sample_id'].map(get_age_group)

    # Save to expected output location
    df.to_csv(OUTPUT_DIR / "cosmic_signatures.csv")
    print(f"Saved COSMIC signatures to {OUTPUT_DIR / 'cosmic_signatures.csv'}")


def generate_signature_figure():
    """Generate stacked bar chart of COSMIC signatures."""
    FIGURES_DIR.mkdir(exist_ok=True)

    df = pd.read_csv(OUTPUT_DIR / "cosmic_signatures.csv", index_col=0)

    # Extract sample IDs for sorting and labeling
    df['_sample_id'] = df.index.map(extract_sample_id)

    # Sort samples: Old first (157, 158, 167), then Young (102, 107, 116)
    # This matches the manuscript figure order
    df = df.sort_values('_sample_id', key=lambda x: x.map(
        lambda sid: (sid in YOUNG_SAMPLES, int(sid) if sid.isdigit() else 0)
    ))

    # Create display labels: ID + age suffix (e.g., "157_O", "102_Y")
    display_labels = []
    for idx in df.index:
        sid = extract_sample_id(idx)
        suffix = '_Y' if sid in YOUNG_SAMPLES else '_O'
        display_labels.append(f"{sid}{suffix}")

    # Drop metadata columns for plotting
    plot_cols = [c for c in df.columns if c not in ['age_group', 'sample_id', '_sample_id']]
    plot_df = df[plot_cols]

    # Plot stacked bar chart
    fig, ax = plt.subplots(figsize=(14, 8))
    plot_df.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')

    ax.set_ylabel('Total Counts', fontsize=18)
    ax.set_xlabel('Samples', fontsize=18)
    ax.set_title('SBS Signature Counts per Sample')
    ax.legend(title='SBS Signatures', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Set x-axis labels with proper ID-based age suffixes
    ax.set_xticklabels(display_labels, fontsize=14, rotation=45, ha='right')
    ax.tick_params(axis='y', labelsize=12)

    plt.tight_layout(rect=[0, 0, 0.85, 1])

    # Save figure
    fig.savefig(FIGURES_DIR / "cosmic_signatures.svg", format='svg', bbox_inches='tight')
    fig.savefig(FIGURES_DIR / "cosmic_signatures.png", format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved COSMIC signatures figure to {FIGURES_DIR}")


def main():
    SIGPROFILER_OUTPUT.mkdir(parents=True, exist_ok=True)

    # Check if SigProfiler already ran
    activities_file = SIGPROFILER_OUTPUT / "Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
    if not activities_file.exists():
        run_sigprofiler()
    else:
        print("SigProfiler output exists, skipping re-run")

    parse_signatures()
    generate_signature_figure()


if __name__ == "__main__":
    main()
