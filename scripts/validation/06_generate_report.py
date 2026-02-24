#!/usr/bin/env python3
"""Generate tiered validation report."""

import json
from datetime import datetime
from pathlib import Path

from config import COMPARISONS_DIR, REPORT_DIR, COMPUTATIONAL_FIGURES, WET_LAB_FIGURES


def load_comparison_results() -> dict:
    """Load all comparison results."""
    results = {
        "numerical": [],
        "visual": [],
        "claims": []
    }

    numerical_path = COMPARISONS_DIR / "numerical" / "numerical_comparison.json"
    if numerical_path.exists():
        with open(numerical_path) as f:
            results["numerical"] = json.load(f)

    visual_path = COMPARISONS_DIR / "visual" / "visual_comparison.json"
    if visual_path.exists():
        with open(visual_path) as f:
            results["visual"] = json.load(f)

    claims_path = COMPARISONS_DIR / "claims" / "claim_verification.json"
    if claims_path.exists():
        with open(claims_path) as f:
            results["claims"] = json.load(f)

    return results


def count_by_status(results: list, status_key: str = "status") -> dict:
    """Count results by status."""
    counts = {"PASS": 0, "REVIEW": 0, "FAIL": 0, "PENDING": 0, "EXACT": 0, "CLOSE": 0, "DIFFERS": 0, "MISSING": 0}
    for r in results:
        status = r.get(status_key) or r.get("recommendation", "PENDING")
        counts[status] = counts.get(status, 0) + 1
    return counts


def determine_overall_status(results: dict) -> str:
    """Determine overall validation status."""
    all_items = results["numerical"] + results["visual"] + results["claims"]

    fail_count = sum(1 for r in all_items if r.get("status") == "FAIL" or r.get("recommendation") == "FAIL")
    review_count = sum(1 for r in all_items if r.get("status") == "REVIEW" or r.get("recommendation") == "REVIEW")

    if fail_count > 0:
        return "FAIL - Requires Investigation"
    elif review_count > 0:
        return "REVIEW - Manual Verification Needed"
    else:
        return "PASS - All Validations Successful"


def generate_report(results: dict) -> str:
    """Generate markdown report."""

    num_counts = count_by_status(results["numerical"])
    vis_counts = count_by_status(results["visual"], "recommendation")
    claim_counts = count_by_status(results["claims"])

    overall_status = determine_overall_status(results)

    report = f"""# Manuscript Validation Report

**Generated:** {datetime.now().strftime("%Y-%m-%d %H:%M")}
**Manuscript:** Revised Manuscript_2-1-26.docx
**Codebase:** ERpBRCA_OlderWomen (refactored)

## Executive Summary

| Tier | Total | PASS/EXACT | REVIEW/CLOSE | FAIL/DIFFERS | PENDING |
|------|-------|------------|--------------|--------------|---------|
| Numerical | {len(results['numerical'])} | {num_counts.get('EXACT', 0)} | {num_counts.get('CLOSE', 0)} | {num_counts.get('DIFFERS', 0)} | {num_counts.get('PENDING', 0) + num_counts.get('MISSING', 0)} |
| Visual | {len(results['visual'])} | {vis_counts.get('PASS', 0)} | {vis_counts.get('REVIEW', 0)} | {vis_counts.get('FAIL', 0)} | {vis_counts.get('PENDING', 0)} |
| Claims | {len(results['claims'])} | {claim_counts.get('PASS', 0)} | {claim_counts.get('REVIEW', 0)} | {claim_counts.get('FAIL', 0)} | {claim_counts.get('PENDING', 0)} |

**Overall Status:** {overall_status}

---

## Tier 1: Numerical Validation

Programmatic comparison of quantitative claims to code outputs.

"""

    if results["numerical"]:
        for r in results["numerical"][:20]:  # Show first 20
            status = r.get("status", "PENDING")
            icon = {"EXACT": "✓", "CLOSE": "~", "DIFFERS": "✗", "PENDING": "?", "MISSING": "?"}.get(status, "?")
            report += f"- {icon} **{r.get('claim_type', 'unknown')}**: {r.get('manuscript_value', 'N/A')} → {r.get('code_value', 'pending')}\n"
        if len(results["numerical"]) > 20:
            report += f"\n*... and {len(results['numerical']) - 20} more claims*\n"
    else:
        report += "*No numerical comparisons completed yet.*\n"

    report += """
---

## Tier 2: Visual Validation

Azure OpenAI GPT-5 vision comparison of figure pairs.

"""

    if results["visual"]:
        for r in results["visual"]:
            status = r.get("recommendation", "PENDING")
            icon = {"PASS": "✓", "REVIEW": "~", "FAIL": "✗"}.get(status, "?")
            report += f"- {icon} **{r.get('figure_id', 'unknown')}**: {status}\n"
            if r.get("differences"):
                for diff in r["differences"][:3]:
                    report += f"  - {diff}\n"
    else:
        report += "*No visual comparisons completed yet.*\n"

    report += """
---

## Tier 3: Claim Verification

Text-to-data verification of manuscript claims.

"""

    if results["claims"]:
        for r in results["claims"][:20]:
            status = r.get("status", "PENDING")
            icon = {"PASS": "✓", "REVIEW": "~", "FAIL": "✗", "PENDING": "?"}.get(status, "?")
            claim_text = r.get("claim", "")[:80]
            report += f"- {icon} {claim_text}...\n"
        if len(results["claims"]) > 20:
            report += f"\n*... and {len(results['claims']) - 20} more claims*\n"
    else:
        report += "*No claim verifications completed yet.*\n"

    report += f"""
---

## Figures Requiring Review

"""

    review_items = [r for r in results["visual"] if r.get("recommendation") == "REVIEW"]
    fail_items = [r for r in results["visual"] if r.get("recommendation") == "FAIL"]

    if fail_items:
        report += "### Failed Validations\n\n"
        for r in fail_items:
            report += f"- **{r.get('figure_id')}**: {r.get('notes', 'See details')}\n"

    if review_items:
        report += "\n### Needs Manual Review\n\n"
        for r in review_items:
            report += f"- **{r.get('figure_id')}**: {r.get('notes', 'See details')}\n"

    if not review_items and not fail_items:
        report += "*No figures flagged for review.*\n"

    report += f"""
---

## Appendix: Scope

### Computational Figures (In Scope)

"""
    for fig_id in COMPUTATIONAL_FIGURES:
        report += f"- {fig_id}\n"

    report += """
### Wet Lab Figures (Out of Scope)

"""
    for fig_id in WET_LAB_FIGURES:
        report += f"- {fig_id}\n"

    return report


def main():
    print("=== Generating Validation Report ===\n")

    results = load_comparison_results()

    print(f"Loaded results:")
    print(f"  Numerical: {len(results['numerical'])}")
    print(f"  Visual: {len(results['visual'])}")
    print(f"  Claims: {len(results['claims'])}")

    report = generate_report(results)

    # Save report
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    report_path = REPORT_DIR / "validation_report.md"

    with open(report_path, "w") as f:
        f.write(report)

    print(f"\nReport saved to: {report_path}")

    # Also save raw results
    results_path = REPORT_DIR / "validation_results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"Raw results: {results_path}")


if __name__ == "__main__":
    main()
