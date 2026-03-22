#!/usr/bin/env python
"""_config.py - Shared configuration for organoid single-cell analysis

Centralizes all external data paths and directory constants.
All scripts import from here to avoid path duplication and drift.

When data locations change (e.g., after GEO upload), update DATA_PATHS here.
"""
from pathlib import Path

# Directory layout (relative to this file)
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / "outputs"
FIGURES_DIR = SCRIPT_DIR / "figures"
CONFIGS_DIR = SCRIPT_DIR / "configs"
LOGS_DIR = SCRIPT_DIR / "logs"

# Create output directories
OUTPUT_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)
LOGS_DIR.mkdir(exist_ok=True)

# External data paths (absolute — these point to data outside the repo)
# Update these when data moves (e.g., after GEO upload)
_NEIL_PROJECT = "/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilOrganoidSingleCell"

DATA_PATHS = {
    # Note: scripts need .../outs/per_sample_outs (Cell Ranger Flex output)
    "pool1_cellranger": f"{_NEIL_PROJECT}/pool1_flex/outs/per_sample_outs",
    "pool2_cellranger": f"{_NEIL_PROJECT}/pool2_flex/outs/per_sample_outs",
    "fastqs_pool1": "/ix1/alee/LO_LAB/General/Lab_Data/20260304_NeilOrganoidSingleCell_Alex/Lee_020326_AL1_FlexPool_1_OES9112A1",
    "fastqs_pool2": "/ix1/alee/LO_LAB/General/Lab_Data/20260304_NeilOrganoidSingleCell_Alex/Lee_020326_AL1_FlexPool_2_OES9112A2",
    "reference": f"{_NEIL_PROJECT}/references/refdata-gex-GRCh38-2024-A-custom",
    "probe_set": f"{_NEIL_PROJECT}/references/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
    "processed_h5ad": f"{_NEIL_PROJECT}/data/processed",
}

# Sample configuration
SAMPLE_IDS = ["OS01", "OS02", "OS03", "OS04", "OS05", "OS06", "OS07"]

SAMPLE_POOL_MAP = {
    "OS01": DATA_PATHS["pool1_cellranger"],
    "OS02": DATA_PATHS["pool1_cellranger"],
    "OS03": DATA_PATHS["pool1_cellranger"],
    "OS04": DATA_PATHS["pool1_cellranger"],
    "OS05": DATA_PATHS["pool2_cellranger"],
    "OS06": DATA_PATHS["pool2_cellranger"],
    "OS07": DATA_PATHS["pool2_cellranger"],
}

# QC thresholds
MIN_GENES = 500
MIN_COUNTS = 1000
MAX_PCT_MITO = 15.0
MIN_CELLS_FOR_SCRUBLET = 50


def h5ad_path(name):
    """Return absolute path to a processed h5ad file."""
    return f"{DATA_PATHS['processed_h5ad']}/{name}"


def check_file_exists(filepath, description="file"):
    """Validate that a file exists, raise informative error if not."""
    p = Path(filepath)
    if not p.exists():
        raise FileNotFoundError(f"ERROR: {description} not found: {filepath}")
    print(f"  Found: {p.name}")
    return p


# Treatment label standardization: ICI → fulv (fulvestrant)
# h5ad files store "E1+ICI"/"E2+ICI" but figures should show "fulv"
TREATMENT_RENAME = {"E1+ICI": "E1+fulv", "E2+ICI": "E2+fulv"}

# Canonical treatment order for all figures (Vehicle first, then E1 group, then E2 group)
TREATMENT_ORDER = [
    "Vehicle",
    "E1", "E1+fulv", "E1+HSD17B7i",
    "E2", "E2+fulv", "E2+HSD17B7i",
]


def standardize_treatments(adata):
    """Rename ICI→fulv in adata.obs['treatment'] for figure clarity."""
    if "treatment" in adata.obs.columns:
        adata.obs["treatment"] = adata.obs["treatment"].replace(TREATMENT_RENAME)
    return adata
