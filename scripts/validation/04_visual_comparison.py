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
