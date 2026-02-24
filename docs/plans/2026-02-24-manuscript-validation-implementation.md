# Manuscript Validation Pipeline Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a pipeline that validates refactored code produces identical scientific conclusions to manuscript figures using tiered comparison (numerical, visual, claims).

**Architecture:** Extract figures from PowerPoint and claims from DOCX, compare against analysis outputs using programmatic checks and Azure OpenAI GPT-5 vision, generate tiered validation report.

**Tech Stack:** Python 3.10+, python-pptx, python-docx, openai (Azure), pandas, numpy

---

## Task 1: Create Validation Directory Structure

**Files:**
- Create: `scripts/validation/__init__.py`
- Create: `scripts/validation/config.py`
- Create: `validation/.gitkeep`

**Step 1: Create directory structure**

```bash
mkdir -p scripts/validation
mkdir -p validation/extracted/figures/main
mkdir -p validation/extracted/figures/supplemental
mkdir -p validation/comparisons/numerical
mkdir -p validation/comparisons/visual
mkdir -p validation/comparisons/claims
mkdir -p validation/report/figures
touch scripts/validation/__init__.py
touch validation/.gitkeep
```

**Step 2: Create config module**

Create `scripts/validation/config.py`:

```python
"""Configuration for manuscript validation pipeline."""

import os
from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
MANUSCRIPT_DIR = PROJECT_ROOT.parent / "Manuscript"
VALIDATION_DIR = PROJECT_ROOT / "validation"
ANALYSIS_DIR = PROJECT_ROOT / "analysis"

# Input files
MAIN_PPTX = MANUSCRIPT_DIR / "Aging_Main Figures_2-1-26.pptx"
SUPP_PPTX = MANUSCRIPT_DIR / "Aging_Supplementary Material_2-1-26.pptx"
MANUSCRIPT_DOCX = MANUSCRIPT_DIR / "Revised Manuscript_2-1-26.docx"

# Output directories
EXTRACTED_DIR = VALIDATION_DIR / "extracted"
FIGURES_DIR = EXTRACTED_DIR / "figures"
COMPARISONS_DIR = VALIDATION_DIR / "comparisons"
REPORT_DIR = VALIDATION_DIR / "report"

# Azure OpenAI configuration
AZURE_OPENAI_ENV_PATH = Path("/ihome/alee/alc376/.azure_openai_api_key")
AZURE_OPENAI_ENDPOINT = "https://alc37-meu8zr6o-eastus2.openai.azure.com/"
AZURE_OPENAI_DEPLOYMENT = "gpt-5-chat"
AZURE_OPENAI_API_VERSION = "2024-12-01-preview"

# Figure scope - computational only
COMPUTATIONAL_FIGURES = {
    "Fig. 2": {"panels": ["A", "B"], "analysis": ["03_rat_wes"], "type": "COMPUTATIONAL"},
    "Fig. 4": {"panels": ["D", "E", "F"], "analysis": ["01_human_bulk_rnaseq"], "type": "MIXED"},
    "Fig. 6": {"panels": ["B", "C", "D"], "analysis": ["04_human_scrnaseq"], "type": "MIXED"},
    "Fig. 7": {"panels": ["A", "B", "C", "D", "E"], "analysis": ["01_human_bulk_rnaseq", "04_human_scrnaseq"], "type": "COMPUTATIONAL"},
    "EDF 2": {"panels": ["A", "B", "C", "D"], "analysis": ["02_rat_snrnaseq", "03_rat_wes", "05_rat_bulk_rnaseq"], "type": "COMPUTATIONAL"},
    "EDF 5": {"panels": ["B"], "analysis": ["01_human_bulk_rnaseq"], "type": "MIXED"},
    "EDF 6": {"panels": ["A", "B"], "analysis": ["01_human_bulk_rnaseq"], "type": "COMPUTATIONAL"},
    "EDF 7": {"panels": ["A", "B", "C", "D"], "analysis": ["01_human_bulk_rnaseq"], "type": "COMPUTATIONAL"},
    "EDF 10": {"panels": ["A", "B"], "analysis": ["04_human_scrnaseq"], "type": "MIXED"},
}

WET_LAB_FIGURES = ["Fig. 1", "Fig. 3", "Fig. 5", "Fig. 8", "EDF 1", "EDF 3", "EDF 4", "EDF 8"]


def load_azure_env():
    """Load Azure OpenAI API key from env file."""
    if AZURE_OPENAI_ENV_PATH.exists():
        with open(AZURE_OPENAI_ENV_PATH) as f:
            for line in f:
                line = line.strip()
                if "=" in line and not line.startswith("#"):
                    key, val = line.split("=", 1)
                    key = key.replace("export ", "").strip()
                    val = val.strip().strip("\"'")
                    os.environ[key] = val
```

**Step 3: Commit**

```bash
git add scripts/validation/ validation/
git commit -m "feat(validation): add directory structure and config"
```

---

## Task 2: Figure Extraction from PowerPoint

**Files:**
- Create: `scripts/validation/01_extract_figures.py`

**Step 1: Create figure extraction script**

Create `scripts/validation/01_extract_figures.py`:

```python
#!/usr/bin/env python3
"""Extract figures from PowerPoint files."""

import argparse
import os
from pathlib import Path
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE

from config import MAIN_PPTX, SUPP_PPTX, FIGURES_DIR


def extract_slide_title(slide) -> str:
    """Extract title or first text from slide."""
    for shape in slide.shapes:
        if shape.has_text_frame:
            text = shape.text_frame.text.strip()
            if text and len(text) < 80:
                return text
    return "untitled"


def sanitize_filename(name: str) -> str:
    """Convert title to safe filename."""
    return "".join(c if c.isalnum() or c in "._- " else "_" for c in name).strip()


def extract_figures_from_pptx(pptx_path: Path, output_dir: Path) -> dict:
    """
    Extract all images from a PowerPoint file.

    Returns:
        dict mapping slide_idx -> list of extracted image paths
    """
    prs = Presentation(pptx_path)
    results = {}

    for slide_idx, slide in enumerate(prs.slides):
        slide_num = slide_idx + 1
        title = extract_slide_title(slide)
        safe_title = sanitize_filename(title)

        slide_dir = output_dir / f"slide_{slide_num:02d}_{safe_title}"
        slide_dir.mkdir(parents=True, exist_ok=True)

        img_paths = []
        img_idx = 0

        for shape in slide.shapes:
            if shape.shape_type == MSO_SHAPE_TYPE.PICTURE:
                image = shape.image
                ext = image.ext
                img_path = slide_dir / f"img_{img_idx:02d}.{ext}"

                with open(img_path, "wb") as f:
                    f.write(image.blob)

                img_paths.append(str(img_path))
                img_idx += 1

        results[slide_num] = {
            "title": title,
            "images": img_paths,
            "image_count": len(img_paths)
        }

        print(f"  Slide {slide_num}: '{title}' - {len(img_paths)} images")

    return results


def main():
    parser = argparse.ArgumentParser(description="Extract figures from PowerPoint")
    parser.add_argument("--main-only", action="store_true", help="Only extract main figures")
    parser.add_argument("--supp-only", action="store_true", help="Only extract supplemental figures")
    args = parser.parse_args()

    results = {"main": {}, "supplemental": {}}

    if not args.supp_only:
        print(f"\n=== Extracting main figures from {MAIN_PPTX.name} ===")
        main_dir = FIGURES_DIR / "main"
        results["main"] = extract_figures_from_pptx(MAIN_PPTX, main_dir)

    if not args.main_only:
        print(f"\n=== Extracting supplemental figures from {SUPP_PPTX.name} ===")
        supp_dir = FIGURES_DIR / "supplemental"
        results["supplemental"] = extract_figures_from_pptx(SUPP_PPTX, supp_dir)

    # Save extraction manifest
    import json
    manifest_path = FIGURES_DIR / "extraction_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\n=== Extraction complete ===")
    print(f"Manifest saved to: {manifest_path}")

    total_main = sum(s["image_count"] for s in results.get("main", {}).values())
    total_supp = sum(s["image_count"] for s in results.get("supplemental", {}).values())
    print(f"Main figures: {total_main} images from {len(results.get('main', {}))} slides")
    print(f"Supplemental figures: {total_supp} images from {len(results.get('supplemental', {}))} slides")


if __name__ == "__main__":
    main()
```

**Step 2: Test extraction**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
python3 scripts/validation/01_extract_figures.py --main-only
```

Expected: Images extracted to `validation/extracted/figures/main/`

**Step 3: Commit**

```bash
git add scripts/validation/01_extract_figures.py
git commit -m "feat(validation): add PowerPoint figure extraction"
```

---

## Task 3: Legend and Claim Extraction from DOCX

**Files:**
- Create: `scripts/validation/02_extract_claims.py`

**Step 1: Create claim extraction script**

Create `scripts/validation/02_extract_claims.py`:

```python
#!/usr/bin/env python3
"""Extract figure legends and quantitative claims from manuscript DOCX."""

import argparse
import json
import re
from pathlib import Path
from docx import Document

from config import MANUSCRIPT_DOCX, EXTRACTED_DIR, COMPUTATIONAL_FIGURES, WET_LAB_FIGURES


def extract_legends(doc: Document) -> dict:
    """
    Extract figure legends from end of manuscript.

    Returns:
        dict mapping figure ID -> legend info
    """
    legends = {}
    in_legends_section = False
    current_fig = None
    current_text = []

    for para in doc.paragraphs:
        text = para.text.strip()

        if "Figures & Figure Legends" in text:
            in_legends_section = True
            continue

        if not in_legends_section:
            continue

        # Check for new figure start
        fig_match = re.match(r"^(Fig\. \d+|Extended Data Fig\.? \d+):", text)
        if fig_match:
            # Save previous figure
            if current_fig and current_text:
                legends[current_fig]["full_text"] = " ".join(current_text)

            # Start new figure
            fig_id = fig_match.group(1)
            # Normalize EDF naming
            fig_id = re.sub(r"Extended Data Fig\.? (\d+)", r"EDF \1", fig_id)

            current_fig = fig_id
            current_text = [text]

            # Determine if computational or wet lab
            fig_type = "COMPUTATIONAL" if fig_id in COMPUTATIONAL_FIGURES else "WET_LAB"

            legends[fig_id] = {
                "title": text.split(":")[1].strip()[:200] if ":" in text else text[:200],
                "type": fig_type,
                "panels": {},
                "full_text": ""
            }
        elif current_fig and text:
            current_text.append(text)

            # Extract panel descriptions (A), (B), etc.
            panel_matches = re.findall(r"\(([A-Z])\)\s*([^(]+?)(?=\([A-Z]\)|$)", text)
            for panel_letter, panel_desc in panel_matches:
                legends[current_fig]["panels"][panel_letter] = panel_desc.strip()[:300]

    # Save last figure
    if current_fig and current_text:
        legends[current_fig]["full_text"] = " ".join(current_text)

    return legends


def extract_quantitative_claims(doc: Document) -> list:
    """
    Extract quantitative claims (p-values, correlations, sample sizes) from text.

    Returns:
        list of claim dicts
    """
    claims = []

    # Patterns to match
    patterns = {
        "p_value": r"p\s*[=<>]\s*[\d.]+(?:e-?\d+)?",
        "fdr": r"FDR\s*[=<>]\s*[\d.]+",
        "correlation": r"(?:Pearson|Spearman)?\s*r\s*=\s*-?[\d.]+",
        "sample_size": r"n\s*=\s*\d+",
        "percentage": r"\d+(?:\.\d+)?%",
        "fold_change": r"\d+(?:\.\d+)?-fold",
    }

    for para_idx, para in enumerate(doc.paragraphs):
        text = para.text.strip()
        if not text:
            continue

        # Find which figure this paragraph references
        fig_refs = re.findall(r"(?:Fig\.|Figure|Extended Data Fig\.?)\s*\d+[A-Z]?", text)

        for pattern_name, pattern in patterns.items():
            matches = re.findall(pattern, text, re.IGNORECASE)
            for match in matches:
                claims.append({
                    "type": pattern_name,
                    "value": match,
                    "context": text[:200],
                    "paragraph_idx": para_idx,
                    "figure_refs": fig_refs
                })

    return claims


def main():
    parser = argparse.ArgumentParser(description="Extract claims from manuscript")
    parser.add_argument("--legends-only", action="store_true")
    parser.add_argument("--claims-only", action="store_true")
    args = parser.parse_args()

    print(f"=== Loading manuscript: {MANUSCRIPT_DOCX.name} ===")
    doc = Document(MANUSCRIPT_DOCX)
    print(f"Total paragraphs: {len(doc.paragraphs)}")

    EXTRACTED_DIR.mkdir(parents=True, exist_ok=True)

    if not args.claims_only:
        print("\n=== Extracting figure legends ===")
        legends = extract_legends(doc)

        legends_path = EXTRACTED_DIR / "legends.json"
        with open(legends_path, "w") as f:
            json.dump(legends, f, indent=2)

        print(f"Extracted {len(legends)} figure legends")
        comp_count = sum(1 for l in legends.values() if l["type"] == "COMPUTATIONAL")
        print(f"  Computational: {comp_count}")
        print(f"  Wet lab: {len(legends) - comp_count}")
        print(f"Saved to: {legends_path}")

    if not args.legends_only:
        print("\n=== Extracting quantitative claims ===")
        claims = extract_quantitative_claims(doc)

        claims_path = EXTRACTED_DIR / "claims.json"
        with open(claims_path, "w") as f:
            json.dump(claims, f, indent=2)

        print(f"Extracted {len(claims)} quantitative claims")
        by_type = {}
        for c in claims:
            by_type[c["type"]] = by_type.get(c["type"], 0) + 1
        for t, count in sorted(by_type.items()):
            print(f"  {t}: {count}")
        print(f"Saved to: {claims_path}")


if __name__ == "__main__":
    main()
```

**Step 2: Test extraction**

```bash
python3 scripts/validation/02_extract_claims.py
```

Expected: `legends.json` and `claims.json` created in `validation/extracted/`

**Step 3: Commit**

```bash
git add scripts/validation/02_extract_claims.py
git commit -m "feat(validation): add DOCX legend and claim extraction"
```

---

## Task 4: Output Status Checker

**Files:**
- Create: `scripts/validation/02b_check_outputs.py`

**Step 1: Create output checker script**

Create `scripts/validation/02b_check_outputs.py`:

```python
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
        "outputs/DE_results/",
    ],
    "03_rat_wes": [
        "outputs/cosmic_signatures.csv",
        "outputs/oncoplot_data.csv",
    ],
    "04_human_scrnaseq": [
        "outputs/seurat_processed.rds",
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
```

**Step 2: Run checker**

```bash
python3 scripts/validation/02b_check_outputs.py
```

**Step 3: Commit**

```bash
git add scripts/validation/02b_check_outputs.py
git commit -m "feat(validation): add output status checker"
```

---

## Task 5: Numerical Comparison

**Files:**
- Create: `scripts/validation/03_numerical_comparison.py`

**Step 1: Create numerical comparison script**

Create `scripts/validation/03_numerical_comparison.py`:

```python
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
```

**Step 2: Run comparison**

```bash
python3 scripts/validation/03_numerical_comparison.py
```

**Step 3: Commit**

```bash
git add scripts/validation/03_numerical_comparison.py
git commit -m "feat(validation): add numerical comparison (tier 1)"
```

---

## Task 6: Visual Comparison with Azure OpenAI

**Files:**
- Create: `scripts/validation/04_visual_comparison.py`

**Step 1: Create visual comparison script**

Create `scripts/validation/04_visual_comparison.py`:

```python
#!/usr/bin/env python3
"""Tier 2: Visual comparison using Azure OpenAI GPT-5 vision."""

import argparse
import base64
import json
import time
from pathlib import Path

from openai import AzureOpenAI

from config import (
    AZURE_OPENAI_ENDPOINT,
    AZURE_OPENAI_DEPLOYMENT,
    AZURE_OPENAI_API_VERSION,
    load_azure_env,
    EXTRACTED_DIR,
    COMPARISONS_DIR,
    COMPUTATIONAL_FIGURES,
)


def encode_image(path: Path) -> str:
    """Encode image to base64."""
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")


def get_image_mime(path: Path) -> str:
    """Get MIME type from file extension."""
    ext = path.suffix.lower()
    return {
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".gif": "image/gif",
    }.get(ext, "image/png")


def compare_figures(
    client: AzureOpenAI,
    manuscript_img: Path,
    regenerated_img: Path,
    figure_id: str,
    panel_desc: str,
) -> dict:
    """
    Compare two figures using GPT-5 vision.

    Returns:
        Comparison result dict
    """
    prompt = f"""You are validating scientific figure reproducibility for {figure_id}.

Figure description from manuscript legend:
{panel_desc}

Compare the two images:
- Image 1: Original manuscript figure
- Image 2: Regenerated from refactored code

Analyze carefully:
1. DATA PATTERNS: Are the same trends, clusters, or distributions visible?
2. NUMERICAL VALUES: Do axis labels, statistics, p-values, sample sizes match?
3. VISUAL ELEMENTS: Are colors, legends, annotations consistent?
4. SCIENTIFIC CONCLUSION: Would the same conclusion be drawn from both?

Extract any visible statistics (p-values, correlations, percentages).

Return ONLY valid JSON (no markdown):
{{
    "match_level": "EXACT" or "SIMILAR" or "DIFFERS",
    "confidence": 0.0 to 1.0,
    "data_patterns_match": true or false,
    "numerical_values_match": true or false,
    "visual_elements_match": true or false,
    "same_conclusion": true or false,
    "differences": ["list", "of", "specific", "differences"],
    "extracted_stats": {{"stat_name": "value"}},
    "recommendation": "PASS" or "REVIEW" or "FAIL",
    "notes": "Brief explanation"
}}"""

    ms_mime = get_image_mime(manuscript_img)
    rg_mime = get_image_mime(regenerated_img)

    response = client.chat.completions.create(
        model=AZURE_OPENAI_DEPLOYMENT,
        messages=[{
            "role": "user",
            "content": [
                {"type": "text", "text": prompt},
                {"type": "image_url", "image_url": {
                    "url": f"data:{ms_mime};base64,{encode_image(manuscript_img)}"
                }},
                {"type": "image_url", "image_url": {
                    "url": f"data:{rg_mime};base64,{encode_image(regenerated_img)}"
                }},
            ]
        }],
        max_tokens=1000,
        temperature=0,
    )

    result_text = response.choices[0].message.content.strip()

    # Parse JSON response
    try:
        # Handle potential markdown code blocks
        if result_text.startswith("```"):
            result_text = result_text.split("```")[1]
            if result_text.startswith("json"):
                result_text = result_text[4:]
        result = json.loads(result_text)
    except json.JSONDecodeError:
        result = {
            "match_level": "ERROR",
            "confidence": 0,
            "raw_response": result_text,
            "recommendation": "REVIEW"
        }

    result["figure_id"] = figure_id
    result["manuscript_image"] = str(manuscript_img)
    result["regenerated_image"] = str(regenerated_img)

    return result


def find_figure_pairs() -> list[tuple]:
    """Find matching manuscript and regenerated figure pairs."""
    # TODO: Implement proper figure matching based on extraction manifest
    # For now, return empty list as placeholder
    pairs = []

    manifest_path = EXTRACTED_DIR / "figures" / "extraction_manifest.json"
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest = json.load(f)
        # Match figures based on slide titles and computational figure list
        # This needs customization based on actual figure organization

    return pairs


def main():
    parser = argparse.ArgumentParser(description="Visual figure comparison")
    parser.add_argument("--dry-run", action="store_true", help="List pairs without calling API")
    parser.add_argument("--figure", type=str, help="Process single figure ID")
    args = parser.parse_args()

    print("=== Tier 2: Visual Comparison (Azure OpenAI GPT-5) ===\n")

    # Load Azure credentials
    load_azure_env()
    import os
    api_key = os.environ.get("AZURE_OPENAI_API_KEY")
    if not api_key:
        print("ERROR: AZURE_OPENAI_API_KEY not found")
        return

    client = AzureOpenAI(
        azure_endpoint=AZURE_OPENAI_ENDPOINT,
        api_version=AZURE_OPENAI_API_VERSION,
        api_key=api_key,
    )

    # Load legends for context
    legends_path = EXTRACTED_DIR / "legends.json"
    legends = {}
    if legends_path.exists():
        with open(legends_path) as f:
            legends = json.load(f)

    # Find figure pairs
    pairs = find_figure_pairs()

    if not pairs:
        print("No figure pairs found. Run extraction scripts first.")
        print("\nTo test with specific images:")
        print("  python3 04_visual_comparison.py --figure 'Fig. 2A'")
        return

    if args.dry_run:
        print(f"Found {len(pairs)} figure pairs:")
        for fig_id, ms_img, rg_img in pairs:
            print(f"  {fig_id}: {ms_img.name} <-> {rg_img.name}")
        return

    # Process pairs
    results = []
    output_dir = COMPARISONS_DIR / "visual"
    output_dir.mkdir(parents=True, exist_ok=True)

    for fig_id, ms_img, rg_img in pairs:
        if args.figure and fig_id != args.figure:
            continue

        print(f"Processing {fig_id}...")

        legend_text = legends.get(fig_id, {}).get("full_text", "No legend available")

        result = compare_figures(client, ms_img, rg_img, fig_id, legend_text[:500])
        results.append(result)

        print(f"  Result: {result.get('recommendation', 'UNKNOWN')}")

        # Rate limiting
        time.sleep(1)

    # Save results
    results_path = output_dir / "visual_comparison.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)

    # Summary
    print(f"\n=== Summary ===")
    for status in ["PASS", "REVIEW", "FAIL"]:
        count = sum(1 for r in results if r.get("recommendation") == status)
        print(f"  {status}: {count}")

    print(f"\nSaved to: {results_path}")


if __name__ == "__main__":
    main()
```

**Step 2: Test with dry run**

```bash
python3 scripts/validation/04_visual_comparison.py --dry-run
```

**Step 3: Commit**

```bash
git add scripts/validation/04_visual_comparison.py
git commit -m "feat(validation): add visual comparison with Azure OpenAI (tier 2)"
```

---

## Task 7: Claim Verification

**Files:**
- Create: `scripts/validation/05_claim_verification.py`

**Step 1: Create claim verification script**

Create `scripts/validation/05_claim_verification.py`:

```python
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
```

**Step 2: Run verification**

```bash
python3 scripts/validation/05_claim_verification.py
```

**Step 3: Commit**

```bash
git add scripts/validation/05_claim_verification.py
git commit -m "feat(validation): add claim verification (tier 3)"
```

---

## Task 8: Report Generator

**Files:**
- Create: `scripts/validation/06_generate_report.py`

**Step 1: Create report generator**

Create `scripts/validation/06_generate_report.py`:

```python
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
```

**Step 2: Run report generation**

```bash
python3 scripts/validation/06_generate_report.py
```

**Step 3: Commit**

```bash
git add scripts/validation/06_generate_report.py
git commit -m "feat(validation): add report generator"
```

---

## Task 9: Master Orchestration Script

**Files:**
- Create: `scripts/validation/run_validation.sh`

**Step 1: Create orchestration script**

Create `scripts/validation/run_validation.sh`:

```bash
#!/bin/bash
# Master script for manuscript validation pipeline

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

cd "$PROJ_DIR"

echo "========================================"
echo "Manuscript Validation Pipeline"
echo "========================================"
echo "Project: $PROJ_DIR"
echo "Started: $(date)"
echo ""

echo "=== Step 1: Extract figures from PowerPoint ==="
python3 scripts/validation/01_extract_figures.py
echo ""

echo "=== Step 2: Extract legends and claims from DOCX ==="
python3 scripts/validation/02_extract_claims.py
echo ""

echo "=== Step 3: Check output status ==="
python3 scripts/validation/02b_check_outputs.py
echo ""

echo "=== Step 4: Numerical comparison (Tier 1) ==="
python3 scripts/validation/03_numerical_comparison.py
echo ""

echo "=== Step 5: Visual comparison (Tier 2) ==="
python3 scripts/validation/04_visual_comparison.py
echo ""

echo "=== Step 6: Claim verification (Tier 3) ==="
python3 scripts/validation/05_claim_verification.py
echo ""

echo "=== Step 7: Generate report ==="
python3 scripts/validation/06_generate_report.py
echo ""

echo "========================================"
echo "Validation Complete"
echo "========================================"
echo "Report: $PROJ_DIR/validation/report/validation_report.md"
echo "Finished: $(date)"
```

**Step 2: Make executable and test**

```bash
chmod +x scripts/validation/run_validation.sh
```

**Step 3: Commit**

```bash
git add scripts/validation/run_validation.sh
git commit -m "feat(validation): add master orchestration script"
```

---

## Task 10: Final Integration Test

**Step 1: Run full pipeline**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/SanghoonCodeReview/ERpBRCA_OlderWomen
./scripts/validation/run_validation.sh
```

**Step 2: Review generated report**

```bash
cat validation/report/validation_report.md
```

**Step 3: Final commit**

```bash
git add validation/
git commit -m "feat(validation): complete manuscript validation pipeline

- Extract figures from PowerPoint (main + supplemental)
- Extract legends and claims from DOCX
- Tier 1: Numerical comparison
- Tier 2: Visual comparison (Azure OpenAI GPT-5)
- Tier 3: Claim verification
- Generate tiered validation report"
```
