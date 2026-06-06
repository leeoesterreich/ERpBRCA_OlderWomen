#!/usr/bin/env python
"""
analysis/07_organoid_single_cell/11_crop_fig5de.py

Render single-subpanel cropped versions of Fig 5D and Fig 5E that match
the legend captions exactly:
  - 5D caption: "UMAP shows dominant clusters following treatment"
  - 5E caption: "Heatmap showing the degree of estrogen response overlap"

The full composite SVGs (04_cluster_composition.svg, 06c_estrogen_overlap_umap.svg)
remain on disk; this script just produces single-panel siblings:
  - figures/heterogeneity/04_cluster_umap_only.{svg,pdf,png}
  - figures/heterogeneity/06c_estrogen_overlap_heatmap_only.{svg,pdf,png}

Replays the exact rendering logic from 10_heterogeneity_analysis.py's
module4 / module6 but only the relevant axes.
"""

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from _config import TREATMENT_COLORS, TREATMENT_ORDER, standardize_treatments

INPUT_H5AD = os.environ.get(
    "MECHANISM_H5AD",
    "data/processed/mechanism_explored.h5ad",
)
FIG_DIR = Path("analysis/07_organoid_single_cell/figures/heterogeneity")
FIG_DIR.mkdir(parents=True, exist_ok=True)


def save(fig, base: Path) -> None:
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")


def crop_5d_umap_only(adata) -> None:
    """Replicates module4_cluster_composition axes[2] (UMAP)."""
    clusters = adata.obs["leiden"].values
    treatments = adata.obs["treatment"].values
    cluster_cats = sorted(adata.obs["leiden"].cat.categories.tolist(), key=lambda x: int(x) if x.isdigit() else x,)

    # Same per-cluster-per-treatment Fisher to get dominant treatment per cluster
    rows = []
    for cl in cluster_cats:
        in_cluster = clusters == cl
        for treat in TREATMENT_ORDER:
            in_treat = treatments == treat
            a = int(np.sum(in_cluster & in_treat))
            b = int(np.sum(in_cluster & ~in_treat))
            c = int(np.sum(~in_cluster & in_treat))
            d = int(np.sum(~in_cluster & ~in_treat))
            odds, _ = stats.fisher_exact([[a, b], [c, d]])
            rows.append({"cluster": cl, "treatment": treat, "odds_ratio": odds})
    enrich_df = pd.DataFrame(rows)

    dominant_treat = enrich_df.loc[enrich_df.groupby("cluster")["odds_ratio"].idxmax()]
    cluster_color_map = {row["cluster"]: TREATMENT_COLORS[row["treatment"]] for _, row in dominant_treat.iterrows()}
    cell_colors = [cluster_color_map.get(c, "#cccccc") for c in clusters]

    fig, ax = plt.subplots(figsize=(7, 6))
    umap = adata.obsm["X_umap"]
    ax.scatter(umap[:, 0], umap[:, 1], c=cell_colors, s=0.5, alpha=0.45, rasterized=True)
    for cl in cluster_cats:
        mask = clusters == cl
        cx, cy = umap[mask, 0].mean(), umap[mask, 1].mean()
        ax.text(
            cx,
            cy,
            cl,
            fontsize=10,
            fontweight="bold",
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7),
        )
    ax.set_title("Clusters Colored by Dominant Treatment", fontsize=11, fontweight="bold")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    save(fig, FIG_DIR / "04_cluster_umap_only")
    plt.close(fig)
    print("Saved:", FIG_DIR / "04_cluster_umap_only.{svg,pdf,png}")


def crop_5e_overlap_heatmap_only(adata) -> None:
    """Replicates module6_estrogen_continuum axes[0] (Pairwise Overlap heatmap)."""
    overlap_csv = Path("analysis/07_organoid_single_cell/outputs/heterogeneity/estrogen_continuum_overlap.csv")
    if overlap_csv.exists():
        overlap_matrix = pd.read_csv(overlap_csv, index_col=0)
        overlap_matrix = overlap_matrix.reindex(index=TREATMENT_ORDER, columns=TREATMENT_ORDER)
    else:
        # Recompute from adata if CSV missing
        er_early = adata.obs["ESTROGEN_RESPONSE_EARLY_score"].values
        er_late = adata.obs["ESTROGEN_RESPONSE_LATE_score"].values
        continuum = (er_early + er_late) / 2
        treatments = adata.obs["treatment"].values
        x_eval = np.linspace(continuum.min() - 0.5, continuum.max() + 0.5, 500)
        kdes = {t: stats.gaussian_kde(continuum[treatments == t], bw_method=0.2)(x_eval) for t in TREATMENT_ORDER}
        overlap_matrix = pd.DataFrame(
            np.zeros((len(TREATMENT_ORDER), len(TREATMENT_ORDER))), index=TREATMENT_ORDER, columns=TREATMENT_ORDER,
        )
        dx = x_eval[1] - x_eval[0]
        for t1 in TREATMENT_ORDER:
            for t2 in TREATMENT_ORDER:
                overlap_matrix.loc[t1, t2] = np.sum(np.minimum(kdes[t1], kdes[t2])) * dx

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(
        overlap_matrix,
        annot=True,
        fmt=".2f",
        cmap="YlOrRd",
        ax=ax,
        square=True,
        vmin=0,
        vmax=1,
        cbar_kws={"label": "KDE Overlap"},
    )
    ax.set_title("Pairwise Estrogen Response Overlap", fontsize=11, fontweight="bold")
    save(fig, FIG_DIR / "06c_estrogen_overlap_heatmap_only")
    plt.close(fig)
    print("Saved:", FIG_DIR / "06c_estrogen_overlap_heatmap_only.{svg,pdf,png}")


def main() -> int:
    print(f"Loading {INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)
    adata = standardize_treatments(adata)
    adata.obs["treatment"] = pd.Categorical(adata.obs["treatment"], categories=TREATMENT_ORDER, ordered=True)
    print(f"  cells={adata.n_obs}  vars={adata.n_vars}")

    crop_5d_umap_only(adata)
    crop_5e_overlap_heatmap_only(adata)
    return 0


if __name__ == "__main__":
    sys.exit(main())
