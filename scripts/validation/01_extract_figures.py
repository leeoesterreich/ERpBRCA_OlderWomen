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
