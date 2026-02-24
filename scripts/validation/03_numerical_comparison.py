#!/usr/bin/env python3
"""Tier 1: Numerical comparison of claims to code outputs."""

import json
import re
import pandas as pd
import numpy as np
from pathlib import Path

from config import EXTRACTED_DIR, ANALYSIS_DIR, COMPARISONS_DIR


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

    for claim in claims:
        manuscript_val = parse_claim_value(claim)

        # TODO: Look up corresponding code output based on figure refs
        # For now, mark as PENDING if we can't find the output
        code_val = None  # Placeholder

        status, diff = compare_values(manuscript_val, code_val, claim["type"])

        results.append({
            "claim_type": claim["type"],
            "manuscript_value": claim["value"],
            "manuscript_numeric": manuscript_val,
            "code_value": code_val,
            "status": status,
            "difference": diff,
            "context": claim["context"][:100],
            "figure_refs": claim.get("figure_refs", [])
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
