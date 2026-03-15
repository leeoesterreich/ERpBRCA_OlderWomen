#!/usr/bin/env python
# analysis/07_organoid_single_cell/09_pathway_enrichment_analysis.py
# GSEA and Enrichr pathway enrichment for E1 vs E1+HSD17B7i
#
# Computes DE from post-probe-fix data, then runs:
# 1. GSEA preranked (Hallmark, KEGG, Reactome)
# 2. Enrichr over-representation (up and down genes separately)
#
# Inputs:
#   - data/processed/mechanism_explored.h5ad (from 08_inhibitor_mechanism_exploration.py)
#
# Outputs:
#   - outputs/inhibitor_mechanism/de_E1_vs_E1_HSD17B7i.csv
#   - outputs/inhibitor_mechanism/gsea_E1_vs_E1_HSD17B7i.csv
#   - outputs/inhibitor_mechanism/enrichr_E1_up_vs_inhibitor.csv
#   - outputs/inhibitor_mechanism/enrichr_E1_down_vs_inhibitor.csv
#   - figures/inhibitor_mechanism/gsea_*.png/.pdf
#   - figures/inhibitor_mechanism/enrichr_*.png/.pdf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc
from datetime import datetime
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False
    print("ERROR: gseapy not installed")

from _config import (
    OUTPUT_DIR, FIGURES_DIR, h5ad_path, check_file_exists
)

# =============================================================================
# CONFIGURATION
# =============================================================================
DATA_PATH = h5ad_path("mechanism_explored.h5ad")
RESULTS_DIR = OUTPUT_DIR / "inhibitor_mechanism"
FIG_DIR = FIGURES_DIR / "inhibitor_mechanism"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 0.25

ENRICHR_LIBRARIES = [
    'MSigDB_Hallmark_2020',
    'KEGG_2021_Human',
    'Reactome_2022',
    'GO_Biological_Process_2023',
    'WikiPathway_2023_Human',
]

GSEA_LIBRARIES = [
    'MSigDB_Hallmark_2020',
    'KEGG_2021_Human',
    'Reactome_2022',
]


def log_msg(msg):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")


# =============================================================================
# STEP 1: Compute DE
# =============================================================================
def compute_de(adata, group1='E1', group2='E1+HSD17B7i'):
    """Run Wilcoxon rank-sum DE between two treatment groups."""
    log_msg(f"Computing DE: {group1} vs {group2}")

    # Subset to the two groups
    mask = adata.obs['treatment'].isin([group1, group2])
    sub = adata[mask].copy()
    log_msg(f"  {group1}: {(sub.obs['treatment'] == group1).sum()} cells")
    log_msg(f"  {group2}: {(sub.obs['treatment'] == group2).sum()} cells")

    # Run DE (group1 as reference, so positive score = higher in group1)
    sc.tl.rank_genes_groups(sub, groupby='treatment', groups=[group1],
                            reference=group2, method='wilcoxon')

    # Extract results
    result = sc.get.rank_genes_groups_df(sub, group=group1)
    result = result.rename(columns={
        'names': 'gene',
        'scores': 'score',
        'pvals': 'pval',
        'pvals_adj': 'padj',
        'logfoldchanges': 'log2FC'
    })

    # Sort by score descending
    result = result.sort_values('score', ascending=False).reset_index(drop=True)

    log_msg(f"  {len(result)} genes tested")
    n_up = ((result['padj'] < PADJ_THRESHOLD) & (result['log2FC'] > LOG2FC_THRESHOLD)).sum()
    n_down = ((result['padj'] < PADJ_THRESHOLD) & (result['log2FC'] < -LOG2FC_THRESHOLD)).sum()
    log_msg(f"  Significant: {n_up} up in {group1}, {n_down} down in {group1}")

    return result


# =============================================================================
# STEP 2: GSEA preranked
# =============================================================================
def run_gsea_preranked(de_df, name):
    """Run GSEA preranked on full DE results."""
    if not GSEAPY_AVAILABLE:
        return None

    log_msg(f"Running GSEA preranked: {name}")

    # Create ranked list using score
    ranked = de_df[['gene', 'score']].dropna()
    ranked = ranked.drop_duplicates(subset='gene')
    ranked = ranked.set_index('gene')['score']
    ranked = ranked.sort_values(ascending=False)

    all_results = []
    for lib in GSEA_LIBRARIES:
        try:
            pre_res = gp.prerank(
                rnk=ranked,
                gene_sets=lib,
                processes=4,
                permutation_num=1000,
                outdir=None,
                no_plot=True,
                seed=42
            )
            if pre_res.res2d is not None and len(pre_res.res2d) > 0:
                results = pre_res.res2d.copy()
                results['Library'] = lib
                all_results.append(results)
                n_sig = (results['FDR q-val'].astype(float) < 0.25).sum()
                log_msg(f"  {lib}: {n_sig} significant (FDR<0.25)")
        except Exception as e:
            log_msg(f"  {lib}: Error - {str(e)[:80]}")

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return None


# =============================================================================
# STEP 3: Enrichr
# =============================================================================
def run_enrichr(gene_list, name):
    """Run Enrichr over-representation analysis."""
    if not GSEAPY_AVAILABLE or len(gene_list) < 5:
        return None

    log_msg(f"Running Enrichr: {name} ({len(gene_list)} genes)")

    all_results = []
    for lib in ENRICHR_LIBRARIES:
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=lib,
                organism='human',
                outdir=None,
                no_plot=True,
                cutoff=0.5
            )
            if enr.results is not None and len(enr.results) > 0:
                results = enr.results.copy()
                results['Library'] = lib
                all_results.append(results)
                n_sig = (results['Adjusted P-value'] < PADJ_THRESHOLD).sum()
                if n_sig > 0:
                    log_msg(f"  {lib}: {n_sig} significant")
        except Exception as e:
            log_msg(f"  {lib}: Error - {str(e)[:80]}")

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    return None


# =============================================================================
# PLOTTING
# =============================================================================
def plot_gsea_dotplot(gsea_df, name, top_n=25):
    """GSEA dot plot: NES on x-axis, terms on y-axis, sized by -log10(FDR)."""
    if gsea_df is None or len(gsea_df) == 0:
        return

    gsea_df = gsea_df.copy()
    gsea_df['FDR q-val'] = gsea_df['FDR q-val'].astype(float)
    gsea_df['NES'] = gsea_df['NES'].astype(float)

    plot_df = gsea_df.nsmallest(top_n, 'FDR q-val').copy()

    # Shorten names
    plot_df['Term_short'] = plot_df['Term'].apply(
        lambda x: x[:55] + '...' if len(str(x)) > 55 else x
    )

    fig, ax = plt.subplots(figsize=(10, max(6, len(plot_df) * 0.4)))

    colors = ['#d32f2f' if nes > 0 else '#1565c0' for nes in plot_df['NES']]
    sizes = -np.log10(plot_df['FDR q-val'].clip(lower=1e-10)) * 25

    ax.scatter(plot_df['NES'], range(len(plot_df)), c=colors, s=sizes, alpha=0.7, edgecolors='k', linewidths=0.5)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df['Term_short'], fontsize=9)
    ax.invert_yaxis()
    ax.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=11)
    ax.set_title(f'GSEA: {name}\n(red = up in E1, blue = down in E1 / up in E1+HSD17B7i)', fontsize=12)

    # Add FDR threshold note
    ax.text(0.02, 0.98, 'Dot size = -log10(FDR)', transform=ax.transAxes,
            fontsize=8, va='top', ha='left', style='italic', color='gray')

    plt.tight_layout()
    fname = f"gsea_{name.lower().replace(' ', '_').replace('+', '_')}"
    plt.savefig(FIG_DIR / f"{fname}.png", dpi=150, bbox_inches='tight')
    plt.savefig(FIG_DIR / f"{fname}.pdf", bbox_inches='tight')
    plt.close()
    log_msg(f"  Saved: {fname}.png/.pdf")


def plot_enrichr_barplot(enr_df, name, top_n=20):
    """Enrichr horizontal bar plot of top enriched terms."""
    if enr_df is None or len(enr_df) == 0:
        return

    sig_df = enr_df[enr_df['Adjusted P-value'] < 0.1].copy()
    if len(sig_df) == 0:
        sig_df = enr_df.head(top_n)

    sig_df = sig_df.nsmallest(top_n, 'Adjusted P-value')
    sig_df['-log10(padj)'] = -np.log10(sig_df['Adjusted P-value'].clip(lower=1e-50))

    sig_df['Term_short'] = sig_df['Term'].apply(
        lambda x: x[:60] + '...' if len(str(x)) > 60 else x
    )

    lib_colors = {
        'MSigDB_Hallmark_2020': '#e74c3c',
        'KEGG_2021_Human': '#3498db',
        'Reactome_2022': '#2ecc71',
        'GO_Biological_Process_2023': '#9b59b6',
        'WikiPathway_2023_Human': '#f39c12'
    }

    fig, ax = plt.subplots(figsize=(12, max(6, len(sig_df) * 0.35)))
    colors = [lib_colors.get(lib, '#333333') for lib in sig_df['Library']]

    ax.barh(range(len(sig_df)), sig_df['-log10(padj)'], color=colors, alpha=0.8)
    ax.set_yticks(range(len(sig_df)))
    ax.set_yticklabels(sig_df['Term_short'], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel('-log10(adjusted p-value)', fontsize=11)
    ax.set_title(f'Enriched Pathways: {name}', fontsize=12)

    ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')

    # Legend for libraries present
    used_libs = sig_df['Library'].unique()
    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=lib_colors.get(l, '#333'), alpha=0.8) for l in used_libs]
    ax.legend(legend_handles, used_libs, loc='lower right', fontsize=8)

    plt.tight_layout()
    fname = f"enrichr_{name.lower().replace(' ', '_').replace('+', '_')}"
    plt.savefig(FIG_DIR / f"{fname}.png", dpi=150, bbox_inches='tight')
    plt.savefig(FIG_DIR / f"{fname}.pdf", bbox_inches='tight')
    plt.close()
    log_msg(f"  Saved: {fname}.png/.pdf")


# =============================================================================
# MAIN
# =============================================================================
def main():
    log_msg("=" * 60)
    log_msg("Pathway Enrichment: E1 vs E1+HSD17B7i")
    log_msg("=" * 60)

    # Load data
    log_msg(f"Loading {DATA_PATH}")
    check_file_exists(DATA_PATH, "Mechanism-explored h5ad")
    adata = sc.read_h5ad(DATA_PATH)
    log_msg(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # --- DE ---
    de_e1 = compute_de(adata, group1='E1', group2='E1+HSD17B7i')
    de_e1.to_csv(RESULTS_DIR / "de_E1_vs_E1_HSD17B7i.csv", index=False)

    # --- GSEA preranked ---
    gsea_e1 = run_gsea_preranked(de_e1, 'E1 vs E1+HSD17B7i')
    if gsea_e1 is not None:
        gsea_e1.to_csv(RESULTS_DIR / "gsea_E1_vs_E1_HSD17B7i.csv", index=False)
        plot_gsea_dotplot(gsea_e1, 'E1 vs E1+HSD17B7i')

    # --- Enrichr: upregulated in E1 (= downregulated by inhibitor) ---
    up_genes = de_e1[(de_e1['padj'] < PADJ_THRESHOLD) & (de_e1['log2FC'] > LOG2FC_THRESHOLD)]['gene'].tolist()
    enr_up = run_enrichr(up_genes, 'E1 Up (Inhibitor Suppresses)')
    if enr_up is not None:
        enr_up.to_csv(RESULTS_DIR / "enrichr_E1_up_vs_inhibitor.csv", index=False)
        plot_enrichr_barplot(enr_up, 'Suppressed by HSD17B7i (E1 context)')

    # --- Enrichr: downregulated in E1 (= upregulated by inhibitor) ---
    down_genes = de_e1[(de_e1['padj'] < PADJ_THRESHOLD) & (de_e1['log2FC'] < -LOG2FC_THRESHOLD)]['gene'].tolist()
    enr_down = run_enrichr(down_genes, 'E1 Down (Inhibitor Activates)')
    if enr_down is not None:
        enr_down.to_csv(RESULTS_DIR / "enrichr_E1_down_vs_inhibitor.csv", index=False)
        plot_enrichr_barplot(enr_down, 'Activated by HSD17B7i (E1 context)')

    # --- Summary ---
    log_msg("\n" + "=" * 60)
    log_msg("COMPLETE")
    log_msg("=" * 60)
    log_msg(f"DE results: {RESULTS_DIR / 'de_E1_vs_E1_HSD17B7i.csv'}")
    if gsea_e1 is not None:
        log_msg(f"GSEA results: {RESULTS_DIR / 'gsea_E1_vs_E1_HSD17B7i.csv'}")
    log_msg(f"Enrichr up: {len(up_genes)} genes queried")
    log_msg(f"Enrichr down: {len(down_genes)} genes queried")
    log_msg(f"Figures: {FIG_DIR}/")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}")
        import traceback
        traceback.print_exc()
        raise
