#!/usr/bin/env python
"""
analysis/06_spatial_biopsies/03_fig7g_immunosuppressive_macrophages.py

Figure 7G - Spatial distribution of Immunosuppressive Macrophages (CD163+)
across n=6 core biopsies from the NCT05914792 trial cohort.

Loads CITEgeist-deconvolved biopsy adata objects and renders one
sc.pl.spatial panel per patient, colored by the
"Immunosuppressive Macrophages (CD163+)" proportion column in adata.obs.

Captures the matplotlib figure manually so both PNG and SVG carry real content.

Output: figures/fig7g/HCC22-088-{P1,P2,P3-S1_A,P4,P5,P6}_immunosuppressive_macrophages_spatial.{svg,png,pdf}
"""

import logging
import os
import pickle
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

BIOPSY_ADATAS_PKL = os.environ.get(
    "BIOPSY_ADATAS_PKL",
    "data/external/biopsy_adatas.pkl",
)
TARGET_COL = "Immunosuppressive Macrophages (CD163+)"
PATIENTS = [
    "HCC22-088-P1-S1",
    "HCC22-088-P2-S1",
    "HCC22-088-P3-S1_A",
    "HCC22-088-P4-S1",
    "HCC22-088-P5-S1",
    "HCC22-088-P6-S1",
]

OUT_DIR = "analysis/06_spatial_biopsies/figures/fig7g"
os.makedirs(OUT_DIR, exist_ok=True)


def main() -> int:
    if not os.path.exists(BIOPSY_ADATAS_PKL):
        log.error("biopsy_adatas.pkl not found: %s", BIOPSY_ADATAS_PKL)
        return 1
    log.info("Loading biopsy adatas from %s", BIOPSY_ADATAS_PKL)
    with open(BIOPSY_ADATAS_PKL, "rb") as f:
        adatas = pickle.load(f)
    log.info("Loaded %d biopsy samples", len(adatas))

    for pid in PATIENTS:
        if pid not in adatas:
            log.warning("Patient %s not in adatas; skipping", pid)
            continue
        adata = adatas[pid]

        # Confirm the target column is present
        if TARGET_COL not in adata.obs.columns:
            available = [c for c in adata.obs.columns if "Immunosuppress" in c or "CD163" in c]
            log.warning("%s: %r missing from obs. Macrophage-ish columns: %s. Skipping.", pid, TARGET_COL, available)
            continue

        log.info("Rendering %s ...", pid)
        # Render with scanpy; do NOT pass save= here (scanpy closes the fig before our savefig)
        sc.pl.spatial(
            adata, color=TARGET_COL, cmap="RdBu_r", vmin=0, vmax=0.6, title=pid, show=False, return_fig=False,
        )
        fig = plt.gcf()
        base = os.path.join(OUT_DIR, f"{pid}_immunosuppressive_macrophages_spatial")
        fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
        fig.savefig(f"{base}.svg", bbox_inches="tight")
        fig.savefig(f"{base}.pdf", bbox_inches="tight")
        plt.close(fig)
        log.info("  -> %s.{png,svg,pdf}", base)

    log.info("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
