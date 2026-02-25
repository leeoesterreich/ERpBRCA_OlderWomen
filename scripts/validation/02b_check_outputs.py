#!/usr/bin/env python3
"""Check which analysis outputs exist for validation."""

import json
from pathlib import Path

from config import ANALYSIS_DIR, EXTRACTED_DIR, COMPUTATIONAL_FIGURES


# Expected outputs per analysis
EXPECTED_OUTPUTS = {
    "01_human_bulk_rnaseq": [
        "outputs/correlation_results.csv",
        "outputs/mica_results.csv",
        "outputs/hsd17b7_correlations.csv",
    ],
    "02_rat_snrnaseq": [
        "outputs/seurat_annotated.rds",
        "../../results/corrected/rat_snrnaseq/DE_results_by_celltype.csv",
    ],
    "03_rat_wes": [
        "outputs/cosmic_signatures.csv",
        "outputs/oncoplot_data.csv",
    ],
    "04_human_scrnaseq": [
        "outputs/seurat_annotated.rds",
        "outputs/tam_analysis.csv",
    ],
    "05_rat_bulk_rnaseq": [
        "outputs/deseq2_results.csv",
        "outputs/pam50_subtypes.csv",
    ],
}


def check_outputs() -> dict:
    """Check which expected outputs exist."""
    results = {
        "found": [],
        "missing": [],
        "by_analysis": {}
    }

    for analysis, expected_files in EXPECTED_OUTPUTS.items():
        analysis_dir = ANALYSIS_DIR / analysis
        analysis_results = {"found": [], "missing": []}

        for expected in expected_files:
            path = analysis_dir / expected
            if path.exists():
                analysis_results["found"].append(expected)
                results["found"].append(f"{analysis}/{expected}")
            else:
                analysis_results["missing"].append(expected)
                results["missing"].append(f"{analysis}/{expected}")

        results["by_analysis"][analysis] = analysis_results

    return results


def check_figures_dir() -> dict:
    """Check which regenerated figures exist."""
    from config import PROJECT_ROOT
    figures_dir = PROJECT_ROOT / "figures" / "by_analysis"

    results = {"found": [], "missing": []}

    for fig_id, fig_info in COMPUTATIONAL_FIGURES.items():
        for analysis in fig_info["analysis"]:
            analysis_fig_dir = figures_dir / analysis.replace("_", "_")
            if analysis_fig_dir.exists():
                pngs = list(analysis_fig_dir.glob("*.png"))
                pdfs = list(analysis_fig_dir.glob("*.pdf"))
                results["found"].extend([str(p) for p in pngs + pdfs])

    return results


def main():
    print("=== Checking Analysis Outputs ===\n")

    output_status = check_outputs()

    print("By Analysis:")
    for analysis, status in output_status["by_analysis"].items():
        found = len(status["found"])
        missing = len(status["missing"])
        total = found + missing
        icon = "✓" if missing == 0 else "✗" if found == 0 else "~"
        print(f"  {icon} {analysis}: {found}/{total} outputs found")
        for m in status["missing"]:
            print(f"      Missing: {m}")

    print(f"\nSummary:")
    print(f"  Found: {len(output_status['found'])}")
    print(f"  Missing: {len(output_status['missing'])}")

    # Save status
    status_path = EXTRACTED_DIR / "output_status.json"
    with open(status_path, "w") as f:
        json.dump(output_status, f, indent=2)
    print(f"\nSaved to: {status_path}")

    # Check figures
    print("\n=== Checking Regenerated Figures ===")
    fig_status = check_figures_dir()
    print(f"  Found: {len(fig_status['found'])} figure files")


if __name__ == "__main__":
    main()
