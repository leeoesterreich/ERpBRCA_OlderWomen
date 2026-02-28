#!/usr/bin/env python
"""
Run CellPhoneDB analysis for Elderly and Young age groups.
Uses CellPhoneDB Python API since CLI may not be in PATH.

Inputs:
    outputs/cellphonedb/counts_*.txt
    outputs/cellphonedb/metadata_*.txt

Outputs:
    outputs/cellphonedb/results_elderly/
    outputs/cellphonedb/results_young/
"""

import os
from pathlib import Path
import pandas as pd

# CellPhoneDB v5 imports
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# Setup paths
SCRIPT_DIR = Path(__file__).parent
INPUT_DIR = SCRIPT_DIR / "outputs" / "cellphonedb"
DB_DIR = INPUT_DIR / "db"
CPDB_FILE = DB_DIR / "cellphonedb.zip"
os.makedirs(INPUT_DIR, exist_ok=True)


def run_cellphonedb(group_name: str):
    """Run CellPhoneDB statistical analysis for a group."""
    print(f"\n{'='*60}")
    print(f"Running CellPhoneDB for {group_name}")
    print(f"{'='*60}")

    counts_file = INPUT_DIR / f"counts_{group_name}.txt"
    meta_file = INPUT_DIR / f"metadata_{group_name}.txt"
    output_dir = INPUT_DIR / f"results_{group_name.lower()}"

    # Check inputs exist
    if not counts_file.exists():
        print(f"ERROR: {counts_file} not found")
        return False
    if not meta_file.exists():
        print(f"ERROR: {meta_file} not found")
        return False

    # Load and validate metadata
    meta_df = pd.read_csv(meta_file, sep='\t')
    print(f"Metadata: {len(meta_df)} cells")
    print(f"Cell types: {meta_df['cell_type'].nunique()}")
    print(f"Cell type counts:\n{meta_df['cell_type'].value_counts()}")

    # Check for empty cell types
    empty_ct = meta_df['cell_type'].isna() | (meta_df['cell_type'] == '')
    if empty_ct.any():
        print(f"WARNING: {empty_ct.sum()} cells have empty cell_type")
        # Filter out empty cell types
        meta_df = meta_df[~empty_ct]
        print(f"Filtered to {len(meta_df)} cells with valid cell types")

    if meta_df['cell_type'].nunique() < 2:
        print("ERROR: Need at least 2 cell types for CellPhoneDB")
        return False

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Run CellPhoneDB statistical analysis
    print(f"\nRunning CellPhoneDB statistical_analysis...")
    print(f"Output: {output_dir}")

    try:
        # Call CellPhoneDB
        cpdb_results = cpdb_statistical_analysis_method.call(
            cpdb_file_path=str(CPDB_FILE),  # Use downloaded database
            meta_file_path=str(meta_file),
            counts_file_path=str(counts_file),
            counts_data='gene_name',
            output_path=str(output_dir),
            threshold=0.1,
            iterations=1000,
            threads=4,
        )

        print(f"\nCellPhoneDB completed successfully!")

        # Check outputs
        expected_files = ['pvalues.csv', 'means.csv', 'significant_means.csv',
                          'deconvoluted.csv']
        for f in expected_files:
            fpath = output_dir / f
            if fpath.exists():
                df = pd.read_csv(fpath)
                print(f"  {f}: {len(df)} rows")
            else:
                print(f"  {f}: NOT FOUND")

        return True

    except Exception as e:
        print(f"ERROR running CellPhoneDB: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    print("=" * 60)
    print("CellPhoneDB Analysis Pipeline")
    print("=" * 60)

    # Run for both age groups
    results = {}
    for group in ["Elderly", "Young"]:
        results[group] = run_cellphonedb(group)

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for group, success in results.items():
        status = "SUCCESS" if success else "FAILED"
        print(f"  {group}: {status}")

    if all(results.values()):
        print("\nAll analyses completed successfully!")
        return 0
    else:
        print("\nSome analyses failed. Check logs above.")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
