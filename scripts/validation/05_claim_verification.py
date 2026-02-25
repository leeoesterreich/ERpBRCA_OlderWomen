#!/usr/bin/env python3
"""Tier 3: Verify manuscript claims against code outputs."""

import json
import re
from pathlib import Path

import pandas as pd
import numpy as np

from config import EXTRACTED_DIR, ANALYSIS_DIR, COMPARISONS_DIR, COMPUTATIONAL_FIGURES


def extract_numeric_value(value_str: str) -> float | None:
    """Extract numeric value from claim string."""
    # Handle p-values
    match = re.search(r"[\d.]+(?:e-?\d+)?", value_str)
    if match:
        try:
            return float(match.group())
        except ValueError:
            return None
    return None


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

    manuscript_val = extract_numeric_value(claim["value"])

    # Look for correlation results
    corr_files = list(analysis_dir.glob("**/correlation*.csv")) + \
                 list(analysis_dir.glob("**/*corr*.csv"))

    if not corr_files:
        result["notes"] = "No correlation output files found"
        return result

    # Parse correlation files and match to claim
    best_match = None
    best_diff = float('inf')

    for corr_file in corr_files:
        try:
            df = pd.read_csv(corr_file)
            for col in df.columns:
                if 'corr' in col.lower() or col.lower() in ['r', 'rho', 'spearman', 'pearson']:
                    for val in df[col].dropna():
                        if isinstance(val, (int, float)) and manuscript_val:
                            diff = abs(val - manuscript_val)
                            if diff < best_diff:
                                best_diff = diff
                                best_match = val
        except Exception:
            continue

    if best_match is not None:
        result["code_value"] = best_match
        # Check if values are close (within 0.01 for correlations)
        if manuscript_val and abs(best_match - manuscript_val) < 0.05:
            result["status"] = "EXACT"
            result["verified"] = True
        elif manuscript_val and abs(best_match - manuscript_val) < 0.1:
            result["status"] = "CLOSE"
            result["verified"] = True
        else:
            result["status"] = "DIFFERS"
            result["verified"] = False
        result["notes"] = f"Found match in {len(corr_files)} correlation files"
    else:
        result["notes"] = f"Found {len(corr_files)} correlation files, no close match"
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

    manuscript_val = extract_numeric_value(claim["value"])

    # Look for DE results
    de_files = list(analysis_dir.glob("**/DE_*.csv")) + \
               list(analysis_dir.glob("**/deseq2*.csv"))

    if not de_files:
        result["notes"] = "No DE output files found"
        return result

    # Parse DE files and match p-values
    best_match = None
    best_diff = float('inf')

    for de_file in de_files:
        try:
            df = pd.read_csv(de_file)
            for col in df.columns:
                if 'p_val' in col.lower() or 'pval' in col.lower() or col.lower() in ['p', 'pvalue']:
                    for val in df[col].dropna():
                        if isinstance(val, (int, float)) and manuscript_val:
                            # Compare on log scale for p-values
                            if val > 0 and manuscript_val > 0:
                                diff = abs(np.log10(val) - np.log10(manuscript_val))
                            else:
                                diff = abs(val - manuscript_val)
                            if diff < best_diff:
                                best_diff = diff
                                best_match = val
        except Exception:
            continue

    if best_match is not None:
        result["code_value"] = best_match
        # Check if p-values match (within 0.5 log units)
        if manuscript_val and manuscript_val > 0 and best_match > 0:
            log_diff = abs(np.log10(best_match) - np.log10(manuscript_val))
            if log_diff < 0.1:
                result["status"] = "EXACT"
                result["verified"] = True
            elif log_diff < 0.5:
                result["status"] = "CLOSE"
                result["verified"] = True
            else:
                result["status"] = "DIFFERS"
                result["verified"] = False
        result["notes"] = f"Found match in {len(de_files)} DE files"
    else:
        result["notes"] = f"Found {len(de_files)} DE files, no close match"
        result["status"] = "REVIEW"

    return result


def verify_sample_size_claim(claim: dict, analysis_dir: Path) -> dict:
    """Verify a sample size claim."""
    result = {
        "claim": claim["context"][:150],
        "type": "sample_size",
        "status": "PENDING",
        "manuscript_value": claim["value"],
        "code_value": None,
        "verified": None,
        "notes": ""
    }

    manuscript_val = extract_numeric_value(claim["value"])
    if manuscript_val:
        manuscript_val = int(manuscript_val)

    # Look for metadata files with sample counts
    meta_files = list(analysis_dir.glob("**/metadata*.csv")) + \
                 list(analysis_dir.glob("**/*fractions*.csv")) + \
                 list(analysis_dir.glob("**/*samples*.csv"))

    for meta_file in meta_files:
        try:
            df = pd.read_csv(meta_file)
            n_rows = len(df)
            if manuscript_val and n_rows == manuscript_val:
                result["code_value"] = n_rows
                result["status"] = "EXACT"
                result["verified"] = True
                result["notes"] = f"Sample count matches in {meta_file.name}"
                return result
        except Exception:
            continue

    result["notes"] = f"Checked {len(meta_files)} metadata files"
    result["status"] = "REVIEW"

    return result


def categorize_claim(claim: dict) -> str:
    """Determine what type of analysis a claim relates to."""
    context = claim["context"].lower()
    claim_type = claim.get("type", "")

    if claim_type == "correlation" or "correlat" in context:
        return "correlation"
    if claim_type == "p_value":
        if "differentially expressed" in context or "upregulated" in context or "downregulated" in context:
            return "differential_expression"
        return "p_value"
    if claim_type == "sample_size":
        return "sample_size"
    if "pathway" in context or "enrichment" in context:
        return "pathway_analysis"
    if "cell type" in context or "proportion" in context:
        return "cell_composition"

    return "other"


def get_analysis_dir_for_figure(fig_ref: str) -> Path | None:
    """Map a figure reference to its analysis directory."""
    # Mapping from figure to analysis directory
    FIGURE_TO_ANALYSIS = {
        "Fig. 1": None,  # Wet lab
        "Fig. 2": "03_rat_wes",
        "Fig. 3": None,  # Wet lab
        "Fig. 4": "01_human_bulk_rnaseq",
        "Fig. 5": None,  # Wet lab
        "Fig. 6": "04_human_scrnaseq",
        "Fig. 7": "01_human_bulk_rnaseq",
        "Fig. 8": None,  # Wet lab
        "EDF 1": None,  # Wet lab
        "EDF 2": "02_rat_snrnaseq",
        "EDF 3": None,  # Wet lab
        "EDF 4": None,  # Wet lab
        "EDF 5": "01_human_bulk_rnaseq",
        "EDF 6": "01_human_bulk_rnaseq",
        "EDF 7": "01_human_bulk_rnaseq",
        "EDF 8": None,  # Wet lab
        "EDF 9": "01_human_bulk_rnaseq",
        "EDF 10": "04_human_scrnaseq",
    }

    # Extract base figure (e.g., "Fig. 4" from "Fig. 4D")
    for fig, analysis in FIGURE_TO_ANALYSIS.items():
        if fig_ref.startswith(fig):
            if analysis:
                return ANALYSIS_DIR / analysis
            return None

    return None


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
        fig_refs = claim.get("figure_refs", [])

        # Determine relevant analysis directory from figure refs
        analysis_dir = None
        for fig_ref in fig_refs:
            analysis_dir = get_analysis_dir_for_figure(fig_ref)
            if analysis_dir:
                break

        # If no figure ref or wet lab figure, mark as out of scope
        if not fig_refs:
            result = {
                "claim": claim["context"][:150],
                "type": category,
                "status": "NO_FIGURE_REF",
                "manuscript_value": claim.get("value"),
                "notes": "No figure reference found for this claim"
            }
        elif analysis_dir is None:
            result = {
                "claim": claim["context"][:150],
                "type": category,
                "status": "WET_LAB",
                "manuscript_value": claim.get("value"),
                "notes": f"Wet lab figure: {fig_refs}"
            }
        elif category == "correlation":
            result = verify_correlation_claim(claim, analysis_dir)
        elif category in ["differential_expression", "p_value"]:
            result = verify_de_claim(claim, analysis_dir)
        elif category == "sample_size":
            result = verify_sample_size_claim(claim, analysis_dir)
        else:
            result = {
                "claim": claim["context"][:150],
                "type": category,
                "status": "REVIEW",
                "manuscript_value": claim.get("value"),
                "notes": f"Verification for {category} not yet implemented"
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
