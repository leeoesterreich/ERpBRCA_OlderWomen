#!/usr/bin/env python3
"""Run SigProfiler to extract COSMIC mutational signatures."""

import os
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Sample age groups
YOUNG_SAMPLES = ['102', '107', '116']
OLD_SAMPLES = ['157', '158', '167']

OUTPUT_DIR = Path(__file__).parent / "outputs"
FIGURES_DIR = Path(__file__).parent / "figures"

# SigProfiler paths
VCF_INPUT = Path("/bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/RatWES_Mutect2_VCF_Input")
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


def parse_signatures():
    """Parse SigProfiler output to CSV format."""
    activities_file = SIGPROFILER_OUTPUT / "Assignment_Solution/Activities/Assignment_Solution_Activities.txt"

    if not activities_file.exists():
        raise FileNotFoundError(f"SigProfiler output not found: {activities_file}")

    df = pd.read_csv(activities_file, sep='\t')
    df = df.set_index('Samples')

    # Drop columns where all values are zero
    df = df.loc[:, (df != 0).any(axis=0)]

    # Add age group annotation
    # Sample IDs 102, 107, 116 are from younger rats per original study design
    young_sample_ids = ['102', '107', '116']
    df['age_group'] = df.index.map(lambda x: 'Young' if x[:3] in young_sample_ids else 'Old')

    # Save to expected output location
    df.to_csv(OUTPUT_DIR / "cosmic_signatures.csv")
    print(f"Saved COSMIC signatures to {OUTPUT_DIR / 'cosmic_signatures.csv'}")


def generate_signature_figure():
    """Generate stacked bar chart of COSMIC signatures."""
    FIGURES_DIR.mkdir(exist_ok=True)

    df = pd.read_csv(OUTPUT_DIR / "cosmic_signatures.csv", index_col=0)

    # Drop age_group column for plotting if present
    if 'age_group' in df.columns:
        df = df.drop(columns=['age_group'])

    # Sort samples: Old first (157, 158, 167), then Young (102, 107, 116)
    sample_order = sorted(df.index, key=lambda x: (x[:3] in YOUNG_SAMPLES, x))
    df = df.reindex(sample_order)

    # Plot stacked bar chart
    fig, ax = plt.subplots(figsize=(14, 8))
    df.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')

    ax.set_ylabel('Total Counts', fontsize=18)
    ax.set_xlabel('Sample', fontsize=18)
    ax.set_title('SBS Signature Counts per Sample')
    ax.legend(title='SBS Signatures', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Relabel x-axis with age group suffix
    new_labels = []
    for label in ax.get_xticklabels():
        sample_id = label.get_text()[:3]
        suffix = '_Y' if sample_id in YOUNG_SAMPLES else '_O'
        new_labels.append(sample_id + suffix)
    ax.set_xticklabels(new_labels, fontsize=14, rotation=45, ha='right')
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
