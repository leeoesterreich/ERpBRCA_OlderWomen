#!/usr/bin/env python3
"""08_compile_report.py — Compile all MICA findings into final report.

Adapted from Jian_MICA/scripts/06_compile_report.py

Inputs:
  - outputs/statistical_audit.json
  - outputs/manuscript_claims.json

Output:
  - outputs/final_report.md
"""
import json
from pathlib import Path
from datetime import date

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "outputs"


def load_json(path):
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return None


def main():
    audit = load_json(OUTPUT_DIR / "statistical_audit.json")
    claims = load_json(OUTPUT_DIR / "manuscript_claims.json")

    # Load Codex figure reviews if present
    codex_reviews = {}
    for f in sorted(OUTPUT_DIR.glob("codex_review_*.txt")):
        fig_name = f.stem.replace("codex_review_", "")
        codex_reviews[fig_name] = f.read_text()

    # Build report
    lines = []
    lines.append(f"# MICA Pipeline Review -- Final Report\n")
    lines.append(f"**Date:** {date.today()}\n")
    lines.append(f"**Project:** Jian Zou MICA aging + breast cancer analysis\n")
    lines.append(f"**Manuscript:** Nature Aging revised submission\n\n")
    lines.append("---\n\n")

    # Section 1: Statistical Audit
    lines.append("## 1. Statistical Audit Findings\n\n")
    if audit:
        for sev in ["critical", "warning", "suggestion"]:
            items = [f for f in audit if f["severity"] == sev]
            if items:
                lines.append(f"### {sev.upper()} ({len(items)})\n\n")
                for item in items:
                    lines.append(f"**{item['category']}: {item['title']}**\n\n")
                    lines.append(f"{item['detail']}\n\n")
    else:
        lines.append("*Audit results not available.*\n\n")

    # Section 2: Manuscript Claims
    lines.append("## 2. Manuscript Claims Extracted\n\n")
    if claims:
        lines.append(f"Total paragraphs: {claims['total_paragraphs']}\n")
        lines.append(f"Quantitative claims: {len(claims['quantitative_claims'])}\n")
        lines.append(f"Directional claims: {len(claims['directional_claims'])}\n")
        lines.append(f"Figures referenced: {', '.join(claims['figures_referenced'])}\n\n")

        lines.append("### Key Quantitative Claims\n\n")
        for c in claims["quantitative_claims"][:20]:
            figs = ", ".join(c["figures"]) if c["figures"] else "none"
            stats = "; ".join(f"{k}: {v}" for k, v in c["statistics"].items())
            lines.append(f"- [{figs}] {stats}\n")
            lines.append(f"  > {c['text'][:200]}...\n\n")
    else:
        lines.append("*Claims extraction not available.*\n\n")

    # Section 3: Figure Reviews
    lines.append("## 3. Figure Visual Reviews (Codex)\n\n")
    if codex_reviews:
        for fig, review in codex_reviews.items():
            lines.append(f"### {fig}\n\n")
            lines.append(f"{review}\n\n")
    else:
        lines.append("*Figure reviews not yet completed.*\n\n")

    # Section 4: Claims Verification Matrix
    lines.append("## 4. Claims Verification Matrix\n\n")
    lines.append("| Claim | Figure | Reproduced? | Statistically Sound? | Biologically Coherent? | Caveated? |\n")
    lines.append("|-------|--------|-------------|---------------------|----------------------|----------|\n")
    lines.append("| EstroGene declines with age | 4H | TBD | WARNING: mean-as-median boxplots | Yes (known biology) | No CI shown |\n")
    lines.append("| HSD17B expression varies by age | 4C | TBD | OK (median heatmap) | Yes | No stats on heatmap |\n")
    lines.append("| Inflammatory pathways increase | 5C | TBD | Check color scale range | Plausible | Curation by expert |\n")
    lines.append("| Immune composition shifts | 5E | TBD | WARNING: mean-as-median, ylim clip | Mixed evidence | No CI shown |\n")
    lines.append("| MICA concordance across cohorts | MICA tables | TBD | CRITICAL: pooled null, no +1 correction | N/A | No FDR reported |\n\n")

    lines.append("## 5. Recommendations\n\n")
    lines.append("*To be filled after all stages complete.*\n")

    report_path = OUTPUT_DIR / "final_report.md"
    report_path.write_text("".join(lines))
    print(f"Final report saved to {report_path}")


if __name__ == "__main__":
    main()
