import argparse
import os
import sys
import pickle
import logging
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
import seaborn as sns
from utils import log_memory_usage
from figure_config import setup_figure_params

# Parse command line arguments
parser = argparse.ArgumentParser(description='Immune pathways analysis')
parser.add_argument('--figures-only', action='store_true',
                    help='Skip computation; regenerate figures from saved results')
args = parser.parse_args()

# Create necessary directories
os.makedirs('logs', exist_ok=True)
os.makedirs('figures', exist_ok=True)
os.makedirs('figures/immune_pathways', exist_ok=True)
os.makedirs('figures/immune_pathways/enrichment', exist_ok=True)
os.makedirs('figures/immune_pathways/prerank', exist_ok=True)
os.makedirs('figures/immune_pathways/dotplots', exist_ok=True)

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/02_immune_pathways.log'),
        logging.StreamHandler(sys.stderr)
    ]
)

# Apply standardized figure parameters (14pt minimum, Arial font)
setup_figure_params()

# Force stdout to flush immediately
sys.stdout.reconfigure(line_buffering=True)

logging.info("Starting immune pathways analysis script")
log_memory_usage()

# Define gene sets for pathway analysis
gene_sets = ['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022']

# Macrophage-related keywords to identify relevant cell types
macrophage_keywords = ['macrophage', 'Macrophage', 'CD163', 'CD14', 'HLA-DR', 'CD11c']

def find_macrophage_layers(adata):
    """Find layers that contain macrophage-related information"""
    macrophage_layers = []
    
    # Check layers
    if hasattr(adata, 'layers') and adata.layers:
        for layer_name in adata.layers.keys():
            if any(keyword in layer_name for keyword in macrophage_keywords):
                macrophage_layers.append(layer_name)
    
    # Check obsm keys
    if hasattr(adata, 'obsm') and adata.obsm:
        for obsm_key in adata.obsm.keys():
            if any(keyword in obsm_key for keyword in macrophage_keywords):
                macrophage_layers.append(obsm_key)
    
    return macrophage_layers

def filter_il_tgf_pathways(results_df):
    """Filter pathways containing IL or TGF in the name"""
    if results_df is None or results_df.empty:
        return results_df
    
    # Create boolean mask for pathways containing IL or TGF
    tgf_mask = results_df['Term'].str.contains('TGF', case=False, na=False)
    interleukin_mask = results_df['Term'].str.contains('interleukin', case=False, na=False)
    
    combined_mask =  tgf_mask | interleukin_mask
    
    filtered_df = results_df[combined_mask].copy()
    
    logging.info(f"Found {len(filtered_df)} pathways containing IL/TGF/interleukin terms out of {len(results_df)} total pathways")
    
    return filtered_df

def run_enrichment_analysis(gene_list, gene_set, sample_name, celltype, direction='up'):
    """Run enrichment analysis and create plots"""
    if len(gene_list) < 10:
        logging.warning(f"Too few genes ({len(gene_list)}) for enrichment analysis in {sample_name} {celltype}")
        return None
    
    try:
        # Run enrichment analysis
        enr = gp.enrichr(gene_list=gene_list,
                        gene_sets=gene_set,
                        outdir=None,  # Don't save intermediate files
                        cutoff=0.05)
        
        if enr.results.empty:
            logging.warning(f"No significant pathways found for {sample_name} {celltype} {gene_set}")
            return None

        # Filter to significant results only (Adjusted P-value < 0.05)
        sig_results = enr.results[enr.results['Adjusted P-value'] < 0.05]

        if sig_results.empty:
            logging.warning(f"No pathways with Adjusted P-value < 0.05 for {sample_name} {celltype} {gene_set}")
            return None

        # Filter for IL/TGF pathways
        filtered_results = filter_il_tgf_pathways(sig_results)
        
        if filtered_results.empty:
            logging.info(f"No IL/TGF pathways found for {sample_name} {celltype} {gene_set}")
            return None
        
        # Create enrichment barplot
        if len(filtered_results) > 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Select top 10 pathways by significance
            top_pathways = filtered_results.head(10)
            
            # Create barplot
            y_pos = np.arange(len(top_pathways))
            bars = ax.barh(y_pos, -np.log10(top_pathways['Adjusted P-value']), 
                          color='steelblue', alpha=0.7)
            
            ax.set_yticks(y_pos)
            ax.set_yticklabels([term[:50] + '...' if len(term) > 50 else term
                               for term in top_pathways['Term']], fontsize=14)
            ax.set_xlabel('-log10(Adjusted P-value)')
            ax.set_title(f'{sample_name} {celltype}\n{gene_set} IL/TGF Pathways ({direction}regulated)')
            
            # Add significance line
            ax.axvline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7)
            
            plt.tight_layout()

            # Save plot as PNG and SVG
            output_path = f"figures/immune_pathways/enrichment/{sample_name}_{celltype}_{gene_set}_{direction}_enrichment.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.savefig(output_path.replace('.png', '.svg'), format='svg', bbox_inches='tight')
            plt.close()

            logging.info(f"Saved enrichment plot: {output_path}")
        
        return filtered_results
        
    except Exception as e:
        logging.error(f"Error in enrichment analysis for {sample_name} {celltype}: {e}")
        return None

def run_prerank_analysis(gene_scores, gene_set, sample_name, celltype):
    """Run GSEA prerank analysis"""
    try:
        # Prepare ranked gene list
        ranked_genes = gene_scores.sort_values(ascending=False)
        
        # Run prerank analysis
        pre_res = gp.prerank(rnk=ranked_genes,
                            gene_sets=gene_set,
                            threads=4,
                            min_size=5,
                            max_size=1000,
                            permutation_num=1000,
                            outdir=None,
                            seed=12345,
                            verbose=False)
        
        if pre_res.res2d.empty:
            logging.warning(f"No significant pathways in prerank for {sample_name} {celltype}")
            return None
        
        # Filter for IL/TGF pathways
        filtered_results = filter_il_tgf_pathways(pre_res.res2d)
        
        if filtered_results.empty:
            logging.info(f"No IL/TGF pathways in prerank for {sample_name} {celltype}")
            return None
        
        # Create prerank plot for top pathway using the prerank object's plot method
        if len(filtered_results) > 0:
            top_pathway = filtered_results.iloc[0]['Term']
            
            try:
                # Use the prerank object's plotting method
                axes = pre_res.plot(terms=[top_pathway], 
                                   figsize=(10, 6),
                                   show_ranking=True)
                
                plt.suptitle(f'{sample_name} {celltype}\nGSEA: {top_pathway}', fontsize=14)

                # Save plot
                output_path = f"figures/immune_pathways/prerank/{sample_name}_{celltype}_{gene_set}_prerank.png"
                plt.savefig(output_path, dpi=300, bbox_inches='tight')
                plt.savefig(output_path.replace('.png', '.svg'), bbox_inches='tight')
                plt.close()
                
                logging.info(f"Saved prerank plot: {output_path}")
                
            except Exception as plot_error:
                logging.warning(f"Could not create prerank plot for {sample_name} {celltype}: {plot_error}")
        
        return filtered_results
        
    except Exception as e:
        logging.error(f"Error in prerank analysis for {sample_name} {celltype}: {e}")
        return None

def create_pathway_dotplot(all_results, output_path):
    """Create dot plot showing pathway activity across samples"""
    if not all_results:
        logging.warning("No results to plot in dotplot")
        return

    # Set Arial font from local file (with fallback)
    import matplotlib.font_manager as fm
    arial_path = '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Arial.ttf'
    if not os.path.exists(arial_path):
        arial_path = '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'
    if os.path.exists(arial_path):
        fm.fontManager.addfont(arial_path)
        font_name = fm.FontProperties(fname=arial_path).get_name()
        plt.rcParams['font.family'] = font_name
    else:
        plt.rcParams['font.family'] = 'sans-serif'

    # Sample name mapping
    sample_map = {
        'HCC22-088-P1-S1': 'S1', 'HCC22-088-P2-S1': 'S2',
        'HCC22-088-P3-S1_A': 'S3', 'HCC22-088-P4-S1': 'S4',
        'HCC22-088-P5-S1': 'S5', 'HCC22-088-P6-S1': 'S6'
    }

    # Combine all results
    combined_df = []
    for (sample, celltype, gene_set), results in all_results.items():
        if results is not None and not results.empty:
            df = results.copy()
            df['Sample'] = sample
            df['Celltype'] = celltype
            df['Gene_set'] = gene_set
            combined_df.append(df)

    if not combined_df:
        logging.warning("No data to combine for dotplot")
        return

    full_df = pd.concat(combined_df, ignore_index=True)

    # Filter to CD163+ macrophages only
    full_df = full_df[full_df['Celltype'].str.contains('CD163', case=False)]

    if full_df.empty:
        logging.warning("No CD163+ macrophage data found")
        return

    # Filter to specific interleukin and TGF-beta pathways
    target_pathways = ['interleukin-4', 'interleukin-10', 'interleukin-13',
                       'interleukin-12', 'interleukin-23', 'interleukin-1', 'interleukin-6',
                       'il-4', 'il-10', 'il-13', 'il-12', 'il-23', 'il-1', 'il-6',
                       'tgf-beta', 'tgf beta', 'transforming growth factor']
    il_mask = full_df['Term'].str.lower().str.contains('|'.join(target_pathways), regex=True)
    full_df = full_df[il_mask]

    if full_df.empty:
        logging.warning("No target interleukin or TGF-beta pathways found")
        return

    # Clean pathway names - remove R-HSA-XXXXXXX patterns
    full_df['Term'] = full_df['Term'].str.replace(r'\s*R-HSA-\d+', '', regex=True)
    
    # Check what columns are available and find p-value column
    logging.info(f"Available columns in results: {list(full_df.columns)}")
    
    # Determine which p-value column to use
    p_value_col = None
    possible_pval_cols = ['Adjusted P-value', 'FDR q-val', 'pval', 'qval', 'Adjusted p-value', 'adj_pval']
    
    for col in possible_pval_cols:
        if col in full_df.columns:
            p_value_col = col
            logging.info(f"Using p-value column: {p_value_col}")
            break
    
    if p_value_col is None:
        logging.error(f"No recognized p-value column found. Available columns: {list(full_df.columns)}")
        return
    
    # Create pivot table for dot plot (only Sample column since we filtered to CD163+ only)
    pivot_df = full_df.pivot_table(
        index='Term',
        columns='Sample',
        values=p_value_col,
        aggfunc='min'
    )

    # Get all unique pathways (no minimum occurrence filter for targeted pathways)
    plot_df = pivot_df.copy()

    if plot_df.empty:
        logging.warning("No pathways to plot")
        return

    # Sort columns by sample order (S1-S6)
    sample_order = ['HCC22-088-P1-S1', 'HCC22-088-P2-S1', 'HCC22-088-P3-S1_A',
                    'HCC22-088-P4-S1', 'HCC22-088-P5-S1', 'HCC22-088-P6-S1']
    available_samples = [s for s in sample_order if s in plot_df.columns]
    plot_df = plot_df[available_samples]

    # Create dot plot with adjusted figure size for full pathway names
    fig, ax = plt.subplots(figsize=(12, max(6, 0.45 * len(plot_df))))

    # Create heatmap
    heatmap = sns.heatmap(-np.log10(plot_df.fillna(1)),
                         cmap='Reds',
                         cbar_kws={'label': '-log10(Adj P-value)'},
                         ax=ax,
                         linewidths=0.5,
                         linecolor='white')

    ax.set_title('Interleukin & TGF-β Pathway Activity in CD163+ Macrophages', fontsize=16)
    ax.set_xlabel('Sample', fontsize=14)
    ax.set_ylabel('Pathway', fontsize=14)

    # Customize x-axis labels using sample_map (S1-S6)
    x_labels = [sample_map.get(col, col) for col in plot_df.columns]
    ax.set_xticklabels(x_labels, rotation=0, ha='center', fontsize=11)

    # Customize y-axis labels - show full pathway names (no truncation)
    ax.set_yticklabels(plot_df.index, rotation=0, fontsize=9)

    plt.subplots_adjust(left=0.42, right=0.97, top=0.90, bottom=0.12)

    # Save as PNG
    plt.savefig(output_path, dpi=600, bbox_inches='tight')

    # Save as SVG
    svg_path = output_path.replace('.png', '.svg')
    plt.savefig(svg_path, dpi=600, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Saved pathway dotplot: {output_path}")


# ============================================================================
# FIGURES-ONLY MODE: Skip computation and regenerate figures from saved results
# ============================================================================
if args.figures_only:
    logging.info("FIGURES-ONLY mode: loading saved results")

    enrichment_pkl = 'figures/immune_pathways/enrichment_results.pkl'
    prerank_pkl = 'figures/immune_pathways/prerank_results.pkl'

    if not os.path.exists(enrichment_pkl) or not os.path.exists(prerank_pkl):
        logging.error("Cannot run figures-only mode: result pickles not found")
        logging.error(f"  Expected: {enrichment_pkl} and {prerank_pkl}")
        logging.error("Run full analysis first to generate intermediate results")
        sys.exit(1)

    with open(enrichment_pkl, 'rb') as f:
        all_enrichment_results = pickle.load(f)
    logging.info(f"Loaded {len(all_enrichment_results)} enrichment results")

    with open(prerank_pkl, 'rb') as f:
        all_prerank_results = pickle.load(f)
    logging.info(f"Loaded {len(all_prerank_results)} prerank results")

    # Create summary dot plots
    logging.info("Creating summary dot plots...")

    if all_enrichment_results:
        create_pathway_dotplot(
            all_enrichment_results,
            'figures/immune_pathways/dotplots/enrichment_summary_dotplot.png'
        )
    else:
        logging.warning("No enrichment results to plot")

    if all_prerank_results:
        create_pathway_dotplot(
            all_prerank_results,
            'figures/immune_pathways/dotplots/prerank_summary_dotplot.png'
        )
    else:
        logging.warning("No prerank results to plot")

    logging.info("Figures-only mode complete!")
    sys.exit(0)


# ============================================================================
# FULL ANALYSIS MODE: Run complete pathway analysis
# ============================================================================

# Load biopsy data (support both plain and gzipped pkl)
import gzip

logging.info("Loading biopsy data...")
pkl_path = os.path.join('data', 'biopsy_adatas.pkl')
pkl_gz_path = pkl_path + '.gz'
if os.path.exists(pkl_gz_path):
    with gzip.open(pkl_gz_path, 'rb') as f:
        biopsy_adatas = pickle.load(f)
elif os.path.exists(pkl_path):
    with open(pkl_path, 'rb') as f:
        biopsy_adatas = pickle.load(f)
else:
    raise FileNotFoundError(f"biopsy_adatas.pkl not found in data/")

logging.info(f"Loaded {len(biopsy_adatas)} biopsy samples")

# Main analysis
logging.info("Starting pathway analysis for macrophages...")

all_enrichment_results = {}
all_prerank_results = {}

for sample_name, adata in biopsy_adatas.items():
    logging.info(f"Processing sample: {sample_name}")
    
    # Find macrophage layers
    macrophage_layers = find_macrophage_layers(adata)
    
    if not macrophage_layers:
        logging.warning(f"No macrophage layers found in {sample_name}")
        continue
    
    logging.info(f"Found macrophage layers: {macrophage_layers}")
    
    # Process each macrophage layer/celltype
    for layer in macrophage_layers:
        celltype_name = layer.replace('_genes_pass1', '').replace('_', ' ')
        
        logging.info(f"Analyzing {celltype_name} in {sample_name}")
        
        # Extract gene expression data for this celltype
        if layer in adata.layers:
            expression_data = adata.layers[layer]
        elif layer in adata.obsm:
            expression_data = adata.obsm[layer]
        else:
            continue
        
        # Convert to DataFrame if needed
        if hasattr(expression_data, 'toarray'):
            expression_data = expression_data.toarray()
        
        # Calculate mean expression per gene
        mean_expression = np.mean(expression_data, axis=0)
        gene_scores = pd.Series(mean_expression, index=adata.var_names)
        
        # Remove zero/negative values and sort
        gene_scores = gene_scores[gene_scores > 0].sort_values(ascending=False)
        
        if len(gene_scores) < 100:
            logging.warning(f"Too few expressed genes ({len(gene_scores)}) in {celltype_name}")
            continue
        
        # Get top upregulated genes for enrichment analysis
        top_genes = gene_scores.head(200).index.tolist()
        
        # Run analysis for each gene set
        for gene_set in gene_sets:
            logging.info(f"Running analysis with {gene_set}")
            
            # Enrichment analysis
            enr_results = run_enrichment_analysis(
                top_genes, gene_set, sample_name, celltype_name, 'up'
            )
            
            if enr_results is not None:
                all_enrichment_results[(sample_name, celltype_name, gene_set)] = enr_results
            
            # Prerank analysis
            prerank_results = run_prerank_analysis(
                gene_scores, gene_set, sample_name, celltype_name
            )
            
            if prerank_results is not None:
                all_prerank_results[(sample_name, celltype_name, gene_set)] = prerank_results

# Create summary dot plots
logging.info("Creating summary dot plots...")

# Create enrichment dotplot
if all_enrichment_results:
    create_pathway_dotplot(
        all_enrichment_results, 
        'figures/immune_pathways/dotplots/enrichment_summary_dotplot.png'
    )
else:
    logging.warning("No enrichment results to plot")

# Create prerank dotplot only if there are results
if all_prerank_results:
    create_pathway_dotplot(
        all_prerank_results,
        'figures/immune_pathways/dotplots/prerank_summary_dotplot.png'
    )
else:
    logging.warning("No prerank results to plot")

# Save results
logging.info("Saving results...")

# Save enrichment results
with open('figures/immune_pathways/enrichment_results.pkl', 'wb') as f:
    pickle.dump(all_enrichment_results, f)

# Save prerank results  
with open('figures/immune_pathways/prerank_results.pkl', 'wb') as f:
    pickle.dump(all_prerank_results, f)

# Create summary report
summary_data = []
for (sample, celltype, gene_set), results in all_enrichment_results.items():
    if results is not None and not results.empty:
        summary_data.append({
            'Sample': sample,
            'Celltype': celltype,
            'Gene_set': gene_set,
            'Num_pathways': len(results),
            'Top_pathway': results.iloc[0]['Term'] if len(results) > 0 else 'None',
            'Top_pvalue': results.iloc[0]['Adjusted P-value'] if len(results) > 0 else 1.0
        })

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv('figures/immune_pathways/analysis_summary.csv', index=False)

logging.info("Analysis complete!")
logging.info(f"Results saved to figures/immune_pathways/")
log_memory_usage()
