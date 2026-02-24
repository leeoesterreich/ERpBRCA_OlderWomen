#!/usr/bin/env python3
"""Tier 3: Verify manuscript claims against code outputs."""

import json
import re
from pathlib import Path

from config import EXTRACTED_DIR, ANALYSIS_DIR, COMPARISONS_DIR


def verify_correlation_claim(claim: dict, analysis_dir: Path) -> dict:
    """Verify a correlation claim against analysis outputs."""
    result = {
        "claim": claim["context"][:150],
        "type": "correlation",
        "status": "PENDING",
        "manuscript_value": claim["value"],
        "code_value": None,
        "verified": None,
        "notes": ""
    }

    # Look for correlation results
    corr_files = list(analysis_dir.glob("**/correlation*.csv")) + \
                 list(analysis_dir.glob("**/*corr*.csv"))

    if not corr_files:
        result["notes"] = "No correlation output files found"
        return result

    # TODO: Parse correlation files and match to claim
    result["notes"] = f"Found {len(corr_files)} correlation files to check"
    result["status"] = "REVIEW"

    return result


def verify_de_claim(claim: dict, analysis_dir: Path) -> dict:
    """Verify a differential expression claim."""
    result = {
        "claim": claim["context"][:150],
        "type": "differential_expression",
        "status": "PENDING",
        "manuscript_value": claim["value"],
        "code_value": None,
        "verified": None,
        "notes": ""
    }

    # Look for DE results
    de_files = list(analysis_dir.glob("**/DE_*.csv")) + \
               list(analysis_dir.glob("**/deseq2*.csv"))

    if not de_files:
        result["notes"] = "No DE output files found"
        return result

    result["notes"] = f"Found {len(de_files)} DE files to check"
    result["status"] = "REVIEW"

    return result


def categorize_claim(claim: dict) -> str:
    """Determine what type of analysis a claim relates to."""
    context = claim["context"].lower()

    if "correlat" in context:
        return "correlation"
    if "differentially expressed" in context or "upregulated" in context or "downregulated" in context:
        return "differential_expression"
    if "pathway" in context or "enrichment" in context:
        return "pathway_analysis"
    if "cell type" in context or "proportion" in context:
        return "cell_composition"

    return "other"


def main():
    print("=== Tier 3: Claim Verification ===\n")

    # Load claims
    claims_path = EXTRACTED_DIR / "claims.json"
    if not claims_path.exists():
        print("ERROR: claims.json not found")
        return

    with open(claims_path) as f:
        claims = json.load(f)

    # Load legends for figure-to-analysis mapping
    legends_path = EXTRACTED_DIR / "legends.json"
    legends = {}
    if legends_path.exists():
        with open(legends_path) as f:
            legends = json.load(f)

    print(f"Loaded {len(claims)} claims")

    results = []

    for claim in claims:
        category = categorize_claim(claim)

        # Determine relevant analysis directory
        fig_refs = claim.get("figure_refs", [])
        # TODO: Map figure refs to analysis directories

        if category == "correlation":
            result = verify_correlation_claim(claim, ANALYSIS_DIR / "01_human_bulk_rnaseq")
        elif category == "differential_expression":
            result = verify_de_claim(claim, ANALYSIS_DIR / "02_rat_snrnaseq")
        else:
            result = {
                "claim": claim["context"][:150],
                "type": category,
                "status": "PENDING",
                "notes": "Verification not yet implemented for this claim type"
            }

        result["figure_refs"] = fig_refs
        results.append(result)

    # Summary
    status_counts = {}
    for r in results:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1

    print("\nResults:")
    for status, count in sorted(status_counts.items()):
        print(f"  {status}: {count}")

    # Save results
    output_dir = COMPARISONS_DIR / "claims"
    output_dir.mkdir(parents=True, exist_ok=True)

    results_path = output_dir / "claim_verification.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nSaved to: {results_path}")


if __name__ == "__main__":
    main()
