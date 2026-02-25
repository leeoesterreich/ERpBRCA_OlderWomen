#!/usr/bin/env python3
"""Tier 1: Numerical comparison of claims to code outputs."""

import json
import re
import pandas as pd
import numpy as np
from pathlib import Path

from config import EXTRACTED_DIR, ANALYSIS_DIR, COMPARISONS_DIR, COMPUTATIONAL_FIGURES


# Mapping of figure references to output files containing verifiable statistics
FIGURE_TO_OUTPUTS = {
    "Fig. 2": {
        "dir": "03_rat_wes",
        "files": ["outputs/cosmic_signatures.csv", "outputs/oncoplot_data.csv"],
    },
    "Fig. 4": {
        "dir": "01_human_bulk_rnaseq",
        "files": ["outputs/correlation_results.csv", "outputs/hsd17b7_correlations.csv"],
    },
    "Fig. 6": {
        "dir": "04_human_scrnaseq",
        "files": ["outputs/cell_fractions.csv", "outputs/tam_analysis.csv"],
    },
    "Fig. 7": {
        "dir": "01_human_bulk_rnaseq",
        "files": ["outputs/mica_results.csv", "outputs/correlation_results.csv"],
    },
    "EDF 2": {
        "dir": "02_rat_snrnaseq",
        "files": ["../../results/corrected/rat_snrnaseq/DE_results_by_celltype.csv"],
    },
    "EDF 5": {
        "dir": "01_human_bulk_rnaseq",
        "files": ["outputs/mica_results.csv"],
    },
    "EDF 6": {
        "dir": "01_human_bulk_rnaseq",
        "files": ["outputs/correlation_results.csv", "outputs/hsd17b7_correlations.csv"],
    },
    "EDF 7": {
        "dir": "01_human_bulk_rnaseq",
        "files": ["outputs/correlation_results.csv"],
    },
}


def load_output_statistics(fig_ref: str) -> dict:
    """Load statistics from output files for a given figure reference."""
    stats = {}

    # Extract base figure (e.g., "Fig. 4" from "Fig. 4D")
    base_fig = None
    for fig in FIGURE_TO_OUTPUTS.keys():
        if fig_ref.startswith(fig):
            base_fig = fig
            break

    if not base_fig:
        return stats

    mapping = FIGURE_TO_OUTPUTS[base_fig]
    analysis_dir = ANALYSIS_DIR / mapping["dir"]

    for file_path in mapping["files"]:
        full_path = analysis_dir / file_path
        if not full_path.exists():
            continue

        try:
            df = pd.read_csv(full_path)

            # Extract p-values
            for col in df.columns:
                if 'p_val' in col.lower() or 'pval' in col.lower() or col.lower() == 'p':
                    for val in df[col].dropna():
                        if isinstance(val, (int, float)):
                            stats[f"p_value_{len(stats)}"] = float(val)

            # Extract correlations
            for col in df.columns:
                if 'corr' in col.lower() or col.lower() == 'r':
                    for val in df[col].dropna():
                        if isinstance(val, (int, float)):
                            stats[f"correlation_{len(stats)}"] = float(val)

            # Extract FDR values
            for col in df.columns:
                if 'fdr' in col.lower() or 'adj' in col.lower():
                    for val in df[col].dropna():
                        if isinstance(val, (int, float)):
                            stats[f"fdr_{len(stats)}"] = float(val)

            # Extract sample counts
            if 'n_cells' in df.columns:
                stats['n_cells_total'] = int(df['n_cells'].sum())

        except Exception as e:
            continue

    return stats


def find_closest_match(manuscript_val: float, stats: dict, claim_type: str) -> tuple:
    """Find the closest matching value in output statistics."""
    if not stats:
        return None, None

    best_match = None
    best_diff = float('inf')
    best_key = None

    for key, val in stats.items():
        if not isinstance(val, (int, float)):
            continue

        # Match p-values with p-values, correlations with correlations, etc.
        if claim_type == "p_value" and "p_value" not in key and "fdr" not in key:
            continue
        if claim_type == "correlation" and "corr" not in key.lower():
            continue

        diff = abs(val - manuscript_val)
        if diff < best_diff:
            best_diff = diff
            best_match = val
            best_key = key

    return best_match, best_key


def parse_claim_value(claim: dict) -> float | None:
    """Extract numeric value from claim."""
    value_str = claim["value"]

    # Handle p-values: p = 0.05, p < 0.001
    if claim["type"] == "p_value":
        match = re.search(r"[\d.]+(?:e-?\d+)?", value_str)
        if match:
            return float(match.group())

    # Handle correlations: r = 0.85
    if claim["type"] == "correlation":
        match = re.search(r"-?[\d.]+", value_str)
        if match:
            return float(match.group())

    # Handle sample sizes: n = 115
    if claim["type"] == "sample_size":
        match = re.search(r"\d+", value_str)
        if match:
            return int(match.group())

    # Handle percentages: 45.3%
    if claim["type"] == "percentage":
        match = re.search(r"[\d.]+", value_str)
        if match:
            return float(match.group())

    return None


def compare_values(manuscript_val: float, code_val: float,
                   claim_type: str) -> tuple[str, float]:
    """
    Compare manuscript claim to code output.

    Returns:
        (status, difference)
        status: EXACT, CLOSE, DIFFERS
    """
    if manuscript_val is None or code_val is None:
        return "MISSING", None

    # For p-values, compare on log scale
    if claim_type == "p_value":
        if manuscript_val == 0 or code_val == 0:
            return "EXACT" if manuscript_val == code_val else "DIFFERS", abs(manuscript_val - code_val)
        log_diff = abs(np.log10(manuscript_val) - np.log10(code_val))
        if log_diff < 0.1:
            return "EXACT", log_diff
        elif log_diff < 0.5:
            return "CLOSE", log_diff
        else:
            return "DIFFERS", log_diff

    # For correlations and percentages
    diff = abs(manuscript_val - code_val)
    rel_diff = diff / max(abs(manuscript_val), 0.001)

    if diff < 0.001 or rel_diff < 0.01:
        return "EXACT", diff
    elif rel_diff < 0.1:
        return "CLOSE", diff
    else:
        return "DIFFERS", diff


def run_numerical_comparison(claims: list) -> list:
    """Run numerical comparison for all claims."""
    results = []

    # Cache for loaded statistics per figure
    stats_cache = {}

    for claim in claims:
        manuscript_val = parse_claim_value(claim)
        fig_refs = claim.get("figure_refs", [])

        # Look up corresponding code output based on figure refs
        code_val = None
        matched_key = None

        for fig_ref in fig_refs:
            # Skip wet lab figures
            if any(wet in fig_ref for wet in ["Fig. 1", "Fig. 3", "Fig. 5", "Fig. 8"]):
                continue

            # Load statistics for this figure (with caching)
            if fig_ref not in stats_cache:
                stats_cache[fig_ref] = load_output_statistics(fig_ref)

            stats = stats_cache[fig_ref]

            if stats and manuscript_val is not None:
                code_val, matched_key = find_closest_match(
                    manuscript_val, stats, claim["type"]
                )
                if code_val is not None:
                    break

        status, diff = compare_values(manuscript_val, code_val, claim["type"])

        # If no figure refs or wet lab figure, mark as OUT_OF_SCOPE
        if not fig_refs or all(any(wet in ref for wet in ["Fig. 1", "Fig. 3", "Fig. 5", "Fig. 8", "EDF 1", "EDF 3", "EDF 4", "EDF 8"]) for ref in fig_refs):
            status = "OUT_OF_SCOPE"

        results.append({
            "claim_type": claim["type"],
            "manuscript_value": claim["value"],
            "manuscript_numeric": manuscript_val,
            "code_value": code_val,
            "matched_key": matched_key,
            "status": status,
            "difference": diff,
            "context": claim["context"][:100],
            "figure_refs": fig_refs
        })

    return results


def main():
    print("=== Tier 1: Numerical Comparison ===\n")

    # Load claims
    claims_path = EXTRACTED_DIR / "claims.json"
    if not claims_path.exists():
        print("ERROR: claims.json not found. Run 02_extract_claims.py first.")
        return

    with open(claims_path) as f:
        claims = json.load(f)

    print(f"Loaded {len(claims)} claims")

    # Run comparison
    results = run_numerical_comparison(claims)

    # Summary
    status_counts = {}
    for r in results:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1

    print("\nResults:")
    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count}")

    # Save results
    output_dir = COMPARISONS_DIR / "numerical"
    output_dir.mkdir(parents=True, exist_ok=True)

    results_path = output_dir / "numerical_comparison.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)

    # Also save as CSV for easy viewing
    df = pd.DataFrame(results)
    csv_path = output_dir / "numerical_comparison.csv"
    df.to_csv(csv_path, index=False)

    print(f"\nSaved to: {results_path}")
    print(f"CSV: {csv_path}")


if __name__ == "__main__":
    main()
