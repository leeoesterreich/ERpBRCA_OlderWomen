#!/usr/bin/env python3
"""Run SigProfiler to extract COSMIC mutational signatures."""

import os
import pandas as pd
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent / "outputs"

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


def main():
    SIGPROFILER_OUTPUT.mkdir(parents=True, exist_ok=True)

    # Check if SigProfiler already ran
    activities_file = SIGPROFILER_OUTPUT / "Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
    if not activities_file.exists():
        run_sigprofiler()
    else:
        print("SigProfiler output exists, skipping re-run")

    parse_signatures()


if __name__ == "__main__":
    main()
