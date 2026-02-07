#!/usr/bin/env python3
"""
Create PowerPoint and Word documents summarizing biostatistical corrections.
"""

import os
from pathlib import Path

# Check for required packages
try:
    from pptx import Presentation
    from pptx.util import Inches, Pt
    from pptx.enum.text import PP_ALIGN
except ImportError:
    print("Installing python-pptx...")
    os.system("pip install python-pptx --quiet")
    from pptx import Presentation
    from pptx.util import Inches, Pt
    from pptx.enum.text import PP_ALIGN

try:
    from docx import Document
    from docx.shared import Inches as DocxInches, Pt as DocxPt
    from docx.enum.text import WD_ALIGN_PARAGRAPH
except ImportError:
    print("Installing python-docx...")
    os.system("pip install python-docx --quiet")
    from docx import Document
    from docx.shared import Inches as DocxInches, Pt as DocxPt
    from docx.enum.text import WD_ALIGN_PARAGRAPH

# Paths
script_dir = Path(__file__).parent
project_root = script_dir.parent
fig_dir = project_root / "figures" / "by_analysis"
output_dir = project_root / "docs"
output_dir.mkdir(exist_ok=True)

print("Creating summary documents...")
print(f"Project root: {project_root}")

# =============================================================================
# Create PowerPoint
# =============================================================================
print("\nCreating PowerPoint presentation...")

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)

def add_title_slide(prs, title, subtitle=""):
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Title
    txBox = slide.shapes.add_textbox(Inches(0.5), Inches(2.5), Inches(12.333), Inches(1.5))
    tf = txBox.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(44)
    p.font.bold = True
    p.alignment = PP_ALIGN.CENTER

    # Subtitle
    if subtitle:
        txBox2 = slide.shapes.add_textbox(Inches(0.5), Inches(4), Inches(12.333), Inches(1))
        tf2 = txBox2.text_frame
        p2 = tf2.paragraphs[0]
        p2.text = subtitle
        p2.font.size = Pt(24)
        p2.alignment = PP_ALIGN.CENTER

    return slide

def add_content_slide(prs, title, bullet_points):
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Title
    txBox = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(12.333), Inches(0.8))
    tf = txBox.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(32)
    p.font.bold = True

    # Content
    txBox2 = slide.shapes.add_textbox(Inches(0.5), Inches(1.2), Inches(12.333), Inches(5.5))
    tf2 = txBox2.text_frame

    for i, point in enumerate(bullet_points):
        if i == 0:
            p = tf2.paragraphs[0]
        else:
            p = tf2.add_paragraph()
        p.text = f"• {point}"
        p.font.size = Pt(20)
        p.space_after = Pt(12)

    return slide

def add_image_slide(prs, title, image_path, caption=""):
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Title
    txBox = slide.shapes.add_textbox(Inches(0.5), Inches(0.2), Inches(12.333), Inches(0.6))
    tf = txBox.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(28)
    p.font.bold = True

    # Image - PDFs need to be converted, so we'll just reference them
    image_path = str(image_path)
    supported_formats = ['.png', '.jpg', '.jpeg', '.gif', '.bmp']

    if os.path.exists(image_path) and any(image_path.lower().endswith(ext) for ext in supported_formats):
        try:
            slide.shapes.add_picture(image_path, Inches(1.5), Inches(0.9), width=Inches(10))
        except Exception as e:
            # Add placeholder text if image fails
            txBox2 = slide.shapes.add_textbox(Inches(1), Inches(2.5), Inches(11), Inches(2))
            tf2 = txBox2.text_frame
            p2 = tf2.paragraphs[0]
            p2.text = f"[See figure: {os.path.basename(image_path)}]"
            p2.font.size = Pt(24)
            p2.font.italic = True
            p2.alignment = PP_ALIGN.CENTER
    else:
        # PDF or missing file - add placeholder
        txBox2 = slide.shapes.add_textbox(Inches(1), Inches(2.5), Inches(11), Inches(2))
        tf2 = txBox2.text_frame
        p2 = tf2.paragraphs[0]
        p2.text = f"[See figure: {os.path.basename(image_path)}]"
        p2.font.size = Pt(24)
        p2.font.italic = True
        p2.alignment = PP_ALIGN.CENTER

        # Add file path info
        p3 = tf2.add_paragraph()
        p3.text = f"Location: figures/by_analysis/{'/'.join(image_path.split('/')[-2:])}"
        p3.font.size = Pt(14)
        p3.alignment = PP_ALIGN.CENTER

    # Caption
    if caption:
        txBox3 = slide.shapes.add_textbox(Inches(0.5), Inches(6.8), Inches(12.333), Inches(0.5))
        tf3 = txBox3.text_frame
        p3 = tf3.paragraphs[0]
        p3.text = caption
        p3.font.size = Pt(14)
        p3.alignment = PP_ALIGN.CENTER

    return slide

def add_table_slide(prs, title, headers, rows):
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Title
    txBox = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(12.333), Inches(0.7))
    tf = txBox.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(28)
    p.font.bold = True

    # Table
    num_rows = len(rows) + 1
    num_cols = len(headers)

    left = Inches(0.5)
    top = Inches(1.2)
    width = Inches(12.333)
    height = Inches(0.5 * num_rows)

    table = slide.shapes.add_table(num_rows, num_cols, left, top, width, height).table

    # Header row
    for i, header in enumerate(headers):
        cell = table.cell(0, i)
        cell.text = header
        cell.text_frame.paragraphs[0].font.bold = True
        cell.text_frame.paragraphs[0].font.size = Pt(14)

    # Data rows
    for row_idx, row in enumerate(rows):
        for col_idx, value in enumerate(row):
            cell = table.cell(row_idx + 1, col_idx)
            cell.text = str(value)
            cell.text_frame.paragraphs[0].font.size = Pt(12)

    return slide

# Slide 1: Title
add_title_slide(prs,
    "ERpBRCA_OlderWomen\nBiostatistical Corrections",
    "Code Review Summary - February 2026")

# Slide 2: Overview
add_content_slide(prs, "Overview of Corrections", [
    "Added FDR (Benjamini-Hochberg) correction for multiple testing",
    "Added missing differential expression analysis for rat snRNA-seq",
    "Added missing differential abundance testing for rat snRNA-seq",
    "Added DoubletFinder for doublet detection",
    "Added random seeds for reproducibility",
    "Fixed Seurat v5 compatibility issues (JoinLayers)"
])

# Slide 3: Human Results Summary
add_table_slide(prs, "Human Bulk RNA-seq: Correlation Analysis",
    ["Metric", "Original", "After FDR"],
    [
        ["Total tests", "50", "50"],
        ["Significant (p < 0.05)", "22", "22"],
        ["Significant (FDR < 0.05)", "N/A", "20"],
        ["Reproducibility", "-", "PASS (100%)"]
    ])

# Slide 4: Human correlation figure - Original
add_image_slide(prs, "Human Bulk RNA-seq: Original Correlations (p < 0.05)",
    fig_dir / "human_bulk_rnaseq" / "correlation_bubbleplot_original.png",
    "22 significant correlations at nominal p < 0.05")

# Slide 5: Human correlation figure - FDR corrected
add_image_slide(prs, "Human Bulk RNA-seq: FDR-Corrected Correlations",
    fig_dir / "human_bulk_rnaseq" / "correlation_bubbleplot_fdr.png",
    "20 significant correlations at FDR < 0.05 (2 lost: HSD17B7 correlations)")

# Slide 6: Rat DE Results
add_table_slide(prs, "Rat snRNA-seq: Differential Expression (NEW)",
    ["Cell Type", "Total DE", "FDR < 0.05", "Up in Aged", "Down in Aged"],
    [
        ["Luminal", "3,482", "3,335", "978", "2,357"],
        ["Basal", "3,802", "2,682", "1,811", "871"],
        ["Fibroblasts", "4,747", "1,631", "1,479", "152"],
        ["Endothelial", "4,308", "1,112", "645", "467"],
        ["TOTAL", "16,339", "8,760", "4,913", "3,847"]
    ])

# Slide 7: Rat volcano plots
add_image_slide(prs, "Rat snRNA-seq: DE Volcano Plots",
    fig_dir / "rat_snrnaseq" / "DE_volcano_combined.png",
    "Differential expression between Young and Aged per cell type")

# Slide 8: Rat DA plots
add_image_slide(prs, "Rat snRNA-seq: Differential Abundance",
    fig_dir / "rat_snrnaseq" / "DA_combined.png",
    "Cell type proportions by age group (n=3 per group)")

# Slide 9: Technical Fixes
add_content_slide(prs, "Technical Fixes Applied", [
    "sys.frame(1)$ofile error: Replaced with commandArgs() for Rscript",
    "future.globals.maxSize: Increased to 4GB for large cell populations",
    "Cluster name 'g' prefix: Fixed Seurat v5 numeric cluster mapping",
    "JoinLayers(): Added before FindAllMarkers for Seurat v5",
    "tryCatch scope: Fixed variable scoping in error handlers"
])

# Slide 10: Conclusions
add_content_slide(prs, "Conclusions", [
    "Core results are REPRODUCIBLE - correlation values identical",
    "91% of human correlations (20/22) remain significant after FDR",
    "Rat snRNA-seq now has complete DE/DA analyses (previously missing)",
    "8,760 FDR-significant DE genes identified in rat mammary tissue",
    "All analyses now include proper multiple testing correction",
    "Random seeds ensure future reproducibility"
])

# Save PowerPoint
pptx_path = output_dir / "biostatistical_corrections_summary.pptx"
prs.save(str(pptx_path))
print(f"PowerPoint saved: {pptx_path}")

# =============================================================================
# Create Word Document
# =============================================================================
print("\nCreating Word document...")

doc = Document()

# Title
title = doc.add_heading("ERpBRCA_OlderWomen: Biostatistical Corrections Report", 0)
title.alignment = WD_ALIGN_PARAGRAPH.CENTER

doc.add_paragraph("Code Review Summary - February 2026")
doc.add_paragraph("")

# Executive Summary
doc.add_heading("Executive Summary", level=1)
doc.add_paragraph(
    "This report documents the biostatistical corrections applied during code review "
    "of the ERpBRCA_OlderWomen analysis pipeline. The corrections ensure proper "
    "multiple testing adjustment, complete statistical analyses, and reproducibility."
)

# Overview of Changes
doc.add_heading("Overview of Changes", level=1)

doc.add_heading("1. FDR Correction", level=2)
doc.add_paragraph(
    "Added Benjamini-Hochberg FDR correction for multiple testing throughout "
    "the analysis pipeline. This is the standard approach for controlling false "
    "discovery rate in genomics studies."
)

doc.add_heading("2. Missing Analyses Added", level=2)
doc.add_paragraph(
    "The original rat snRNA-seq analysis was missing critical statistical tests:"
)
bullets = doc.add_paragraph()
bullets.add_run("• Differential expression analysis (Young vs Aged per cell type)\n")
bullets.add_run("• Differential abundance testing (cell type proportions)\n")
bullets.add_run("• DoubletFinder for quality control")

doc.add_heading("3. Reproducibility Improvements", level=2)
doc.add_paragraph(
    "Random seeds (12345) were added to all analyses to ensure reproducibility. "
    "Script path detection was fixed for Rscript execution compatibility."
)

# Human Bulk RNA-seq Results
doc.add_heading("Human Bulk RNA-seq Results", level=1)

doc.add_heading("Correlation Analysis", level=2)
table = doc.add_table(rows=5, cols=3)
table.style = 'Table Grid'
headers = table.rows[0].cells
headers[0].text = "Metric"
headers[1].text = "Original"
headers[2].text = "After FDR"

data = [
    ("Total tests", "50", "50"),
    ("Significant (p < 0.05)", "22", "22"),
    ("Significant (FDR < 0.05)", "N/A", "20"),
    ("Reproducibility", "-", "PASS")
]

for i, (metric, orig, fdr) in enumerate(data):
    row = table.rows[i + 1].cells
    row[0].text = metric
    row[1].text = orig
    row[2].text = fdr

doc.add_paragraph("")
doc.add_paragraph(
    "Key finding: All 22 nominally significant correlations are identical between "
    "original and corrected analyses, confirming reproducibility. After FDR correction, "
    "2 correlations (HSD17B7 × Estrogen pathways) lost significance due to borderline "
    "p-values (0.026, 0.036)."
)

# Rat snRNA-seq Results
doc.add_heading("Rat snRNA-seq Results", level=1)

doc.add_heading("Differential Expression Analysis (NEW)", level=2)
doc.add_paragraph(
    "This analysis was completely missing from the original pipeline. We now have "
    "comprehensive DE testing between Young and Aged samples for each cell type."
)

table2 = doc.add_table(rows=6, cols=5)
table2.style = 'Table Grid'
headers2 = table2.rows[0].cells
headers2[0].text = "Cell Type"
headers2[1].text = "Total DE"
headers2[2].text = "FDR < 0.05"
headers2[3].text = "Up in Aged"
headers2[4].text = "Down in Aged"

de_data = [
    ("Luminal", "3,482", "3,335", "978", "2,357"),
    ("Basal", "3,802", "2,682", "1,811", "871"),
    ("Fibroblasts", "4,747", "1,631", "1,479", "152"),
    ("Endothelial", "4,308", "1,112", "645", "467"),
    ("TOTAL", "16,339", "8,760", "4,913", "3,847")
]

for i, row_data in enumerate(de_data):
    row = table2.rows[i + 1].cells
    for j, val in enumerate(row_data):
        row[j].text = val

doc.add_paragraph("")

doc.add_heading("Biological Insights", level=2)
bullets2 = doc.add_paragraph()
bullets2.add_run("• Luminal cells show predominantly downregulation in aging (2,357 down vs 978 up)\n")
bullets2.add_run("• Fibroblasts show predominantly upregulation in aging (1,479 up vs 152 down)\n")
bullets2.add_run("• Basal cells show mixed regulation with more upregulation\n")
bullets2.add_run("• Total of 8,760 genes show significant age-related changes")

doc.add_heading("Differential Abundance Analysis (NEW)", level=2)
doc.add_paragraph(
    "Propeller testing was performed to compare cell type proportions between age groups. "
    "No significant differences were detected at FDR < 0.05, though this is expected given "
    "the limited sample size (n=3 per group)."
)

# Technical Fixes
doc.add_heading("Technical Fixes", level=1)

fixes = [
    ("sys.frame(1)$ofile error",
     "Replaced with commandArgs() function for proper script path detection when run via Rscript"),
    ("future.globals.maxSize",
     "Increased from 500MB to 4GB to handle large cell populations in Seurat"),
    ("Cluster name 'g' prefix",
     "Fixed Seurat v5 behavior that prepends 'g' to numeric cluster names"),
    ("JoinLayers()",
     "Added before FindAllMarkers and FindMarkers for Seurat v5 compatibility"),
    ("tryCatch scope",
     "Fixed variable scoping so error handlers properly set outer variables")
]

for fix, desc in fixes:
    p = doc.add_paragraph()
    p.add_run(f"{fix}: ").bold = True
    p.add_run(desc)

# Conclusions
doc.add_heading("Conclusions", level=1)

doc.add_paragraph(
    "The biostatistical corrections successfully address all issues identified "
    "during code review:"
)

conclusions = doc.add_paragraph()
conclusions.add_run("1. ").bold = True
conclusions.add_run("Reproducibility confirmed - correlation values are identical\n")
conclusions.add_run("2. ").bold = True
conclusions.add_run("91% of human correlations remain significant after FDR correction\n")
conclusions.add_run("3. ").bold = True
conclusions.add_run("Rat snRNA-seq now has complete DE/DA analyses\n")
conclusions.add_run("4. ").bold = True
conclusions.add_run("8,760 FDR-significant DE genes identified across 4 cell types\n")
conclusions.add_run("5. ").bold = True
conclusions.add_run("All technical issues resolved for Seurat v5 compatibility")

doc.add_paragraph("")
doc.add_paragraph(
    "The core biological findings are robust and the analysis pipeline now meets "
    "standards for biostatistical rigor."
)

# Save Word document
docx_path = output_dir / "biostatistical_corrections_report.docx"
doc.save(str(docx_path))
print(f"Word document saved: {docx_path}")

print("\n=== Documents created successfully ===")
print(f"PowerPoint: {pptx_path}")
print(f"Word: {docx_path}")
