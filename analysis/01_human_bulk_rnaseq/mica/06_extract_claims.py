#!/usr/bin/env python3
"""06_extract_claims.py — Extract quantitative and directional claims from manuscript.

Parses the Nature Aging manuscript .docx and extracts:
- Quantitative claims (p-values, effect sizes, fold changes, sample sizes)
- Directional claims (increases/decreases with age, concordant across cohorts)
- Maps each claim to its supporting figure/table

Output: outputs/manuscript_claims.json

Adapted from Jian_MICA/scripts/04_extract_claims.py
"""
import json
import re
import sys
from pathlib import Path

try:
    from docx import Document
except ImportError:
    sys.exit("ERROR: python-docx not installed. Run: pip install python-docx")

PROJ_DIR = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview")
MANUSCRIPT = PROJ_DIR / "Manuscript" / "Revised Manuscript_2-1-26.docx"

# Output goes to monorepo mica/outputs/
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT = SCRIPT_DIR / "outputs" / "manuscript_claims.json"


def extract_text(docx_path):
    """Extract all paragraph text from a .docx file."""
    doc = Document(str(docx_path))
    return [p.text.strip() for p in doc.paragraphs if p.text.strip()]


def find_quantitative_claims(paragraphs):
    """Find paragraphs containing p-values, effect sizes, or statistical measures."""
    claims = []
    patterns = {
        "p_value": r'[Pp]\s*[<>=]\s*[\d.]+(?:\s*[x]\s*10\s*[-]\s*\d+)?',
        "fold_change": r'\d+\.?\d*\s*[-]?\s*fold',
        "effect_size": r"(?:Cohen'?s?\s*d|Hedge'?s?\s*g|r\s*=)\s*[\d.]+",
        "sample_size": r'[Nn]\s*=\s*\d+',
        "fdr": r'(?:FDR|q)\s*[<>=]\s*[\d.]+',
        "confidence_interval": r'\d+%?\s*(?:CI|confidence interval)',
    }

    for i, para in enumerate(paragraphs):
        found = {}
        for name, pattern in patterns.items():
            matches = re.findall(pattern, para)
            if matches:
                found[name] = matches
        if found:
            fig_match = re.findall(r'(?:Fig(?:ure)?\.?\s*\d+[A-Za-z]?)', para)
            claims.append({
                "paragraph_index": i,
                "text": para[:300],
                "statistics": found,
                "figures": fig_match if fig_match else [],
            })
    return claims


def find_directional_claims(paragraphs):
    """Find paragraphs with directional biological claims."""
    claims = []
    direction_patterns = [
        r'(?:increas|decreas|declin|elevat|reduc|diminish|enrich|deplet)\w*\s+(?:with|in|among)\s+(?:age|aging|elderly|older)',
        r'(?:age|aging)[\w\s]*(?:associat|correlat|relat)\w*\s+(?:with)',
        r'concordan\w+\s+(?:across|between|among)',
        r'(?:young|middle|elderly)[\w\s]*(?:higher|lower|greater|less|more|fewer)',
    ]

    for i, para in enumerate(paragraphs):
        matches = []
        for pattern in direction_patterns:
            found = re.findall(pattern, para, re.IGNORECASE)
            if found:
                matches.extend(found)
        if matches:
            fig_match = re.findall(r'(?:Fig(?:ure)?\.?\s*\d+[A-Za-z]?)', para)
            claims.append({
                "paragraph_index": i,
                "text": para[:300],
                "directional_phrases": matches,
                "figures": fig_match if fig_match else [],
            })
    return claims


def main():
    if not MANUSCRIPT.exists():
        sys.exit(f"ERROR: Manuscript not found at {MANUSCRIPT}")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    print(f"Parsing: {MANUSCRIPT}")
    paragraphs = extract_text(MANUSCRIPT)
    print(f"Total paragraphs: {len(paragraphs)}")

    quant_claims = find_quantitative_claims(paragraphs)
    dir_claims = find_directional_claims(paragraphs)

    print(f"Quantitative claims found: {len(quant_claims)}")
    print(f"Directional claims found: {len(dir_claims)}")

    all_figs = set()
    for c in quant_claims + dir_claims:
        all_figs.update(c.get("figures", []))
    print(f"Figures referenced: {sorted(all_figs)}")

    result = {
        "manuscript": str(MANUSCRIPT),
        "total_paragraphs": len(paragraphs),
        "quantitative_claims": quant_claims,
        "directional_claims": dir_claims,
        "figures_referenced": sorted(all_figs),
    }

    with open(OUTPUT, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved to {OUTPUT}")


if __name__ == "__main__":
    main()
