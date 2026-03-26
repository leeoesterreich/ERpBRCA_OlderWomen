#!/usr/bin/env python3
"""
Run CellPhoneDB v4 analysis for Elderly and Young age groups.
Uses CellPhoneDB v4 API (different from v5).

Inputs:
    outputs/cellphonedb/counts_*.txt
    outputs/cellphonedb/metadata_*.txt

Outputs:
    outputs/cellphonedb_v4/results_elderly/
    outputs/cellphonedb_v4/results_young/
"""

import os
import sys
from pathlib import Path

# Get script directory
SCRIPT_DIR = Path(__file__).parent.resolve()
INPUT_DIR = SCRIPT_DIR / "outputs" / "cellphonedb"
OUTPUT_DIR = SCRIPT_DIR / "outputs" / "cellphonedb_v4"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# CellPhoneDB v4 database
CPDB_FILE = OUTPUT_DIR / "cellphonedb.zip"

def run_cellphonedb_v4(group_name: str):
    """Run CellPhoneDB v4 statistical analysis for a group."""

    counts_file = INPUT_DIR / f"counts_{group_name}.txt"
    meta_file = INPUT_DIR / f"metadata_{group_name}.txt"
    output_path = OUTPUT_DIR / f"results_{group_name.lower()}"

    print(f"\nRunning CellPhoneDB v4 for {group_name}")
    print(f"  Counts: {counts_file}")
    print(f"  Metadata: {meta_file}")
    print(f"  Output: {output_path}")

    if not counts_file.exists():
        print(f"ERROR: Counts file not found: {counts_file}")
        return False

    if not meta_file.exists():
        print(f"ERROR: Metadata file not found: {meta_file}")
        return False

    # CellPhoneDB v4 uses the cellphonedb command line or the method API
    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)

    try:
        # Run CellPhoneDB v4 statistical analysis
        # Note: v4 API is slightly different from v5
        print(f"\nRunning CellPhoneDB v4 statistical_analysis...")

        cpdb_statistical_analysis_method.call(
            cpdb_file_path=str(CPDB_FILE),
            meta_file_path=str(meta_file),
            counts_file_path=str(counts_file),
            counts_data='hgnc_symbol',
            output_path=str(output_path),
            threshold=0.1,
            threads=4,
            debug_seed=42,
            result_precision=3,
        )

        print(f"\nCellPhoneDB v4 completed for {group_name}!")
        return True

    except Exception as e:
        print(f"ERROR running CellPhoneDB v4: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    print("=" * 60)
    print("CellPhoneDB v4 Analysis Pipeline")
    print("=" * 60)

    # Check cellphonedb version
    try:
        import pkg_resources
        version = pkg_resources.get_distribution("cellphonedb").version
        print(f"CellPhoneDB version: {version}")
    except Exception:
        print("CellPhoneDB version: 4.x (version check failed)")

    results = {}
    for group in ["Elderly", "Young"]:
        results[group] = run_cellphonedb_v4(group)

    print("\n" + "=" * 60)
    print("Summary:")
    for group, success in results.items():
        status = "SUCCESS" if success else "FAILED"
        print(f"  {group}: {status}")
    print("=" * 60)

    # Return non-zero if any failed
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
