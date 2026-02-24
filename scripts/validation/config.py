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
