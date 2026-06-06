# PROVENANCE: Ported 2026-06-05 from
#   /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/utils.py
# Provides plot_single_group_pathway (line ~1860) which renders Fig 7I-style
# sender top-pathways bars. Full CITEgeist utils retained for reproducibility.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import os
import logging
import psutil
from scipy import stats


def volcano_plot(adata, comparator):
    # Extract the results from rank genes groups analysis
    results = adata.uns['rank_genes_groups']
    groups = results['names'].dtype.names

    # Create dataframe with DEG results
    deg_df = pd.DataFrame(
        {key: results[key][groups[0]]
        for key in ['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']}
    )

    # Calculate -log10 of adjusted p-values
    deg_df['-log10(pvals_adj)'] = -np.log10(deg_df['pvals_adj'])

    # Set font sizes for all elements
    plt.rcParams.update({'font.size': 18})
    
    # Create volcano plot
    plt.figure(figsize=(12, 8))

    # Plot all points in gray
    plt.scatter(deg_df['logfoldchanges'], deg_df['-log10(pvals_adj)'], 
            c='gray', alpha=0.5, s=20)

    # Get significant genes (padj < 0.05)
    sig_genes = deg_df[deg_df['pvals_adj'] < 0.05].copy()

    # Split into upregulated (positive logFC > 1) and downregulated (negative logFC < -1)
    up_genes = sig_genes[sig_genes['logfoldchanges'] > 1]  
    down_genes = sig_genes[sig_genes['logfoldchanges'] < -1]  

    # Plot upregulated in red and downregulated in blue
    plt.scatter(up_genes['logfoldchanges'], up_genes['-log10(pvals_adj)'],
            c='red', alpha=0.7, s=20, label='Upregulated')
    plt.scatter(down_genes['logfoldchanges'], down_genes['-log10(pvals_adj)'], 
            c='blue', alpha=0.7, s=20, label='Downregulated')

    # Find genes to label - separate tracking for up and down regulated genes
    up_texts = []
    down_texts = []
    
    # For upregulated genes
    if not up_genes.empty:
        # Get top 10 by fold change
        top_fc_up = up_genes.nlargest(10, 'logfoldchanges')
        # Get top 10 by p-value
        top_pval_up = up_genes.nlargest(10, '-log10(pvals_adj)')
        # Find intersection of gene names
        up_genes_to_label = set(top_fc_up['names']).intersection(set(top_pval_up['names']))
        # If intersection is empty, take top 5 from each
        if not up_genes_to_label:
            up_genes_to_label = set(top_fc_up['names'].head(5)).union(set(top_pval_up['names'].head(5)))
        
        # Label the selected genes
        for gene_name in up_genes_to_label:
            gene = up_genes[up_genes['names'] == gene_name].iloc[0]
            up_texts.append(plt.text(gene['logfoldchanges'], 
                                gene['-log10(pvals_adj)'],
                                gene['names'],
                                fontsize=18))
    
    # For downregulated genes
    if not down_genes.empty:
        # Get top 10 by absolute fold change
        down_genes['abs_logfc'] = abs(down_genes['logfoldchanges'])
        top_fc_down = down_genes.nlargest(10, 'abs_logfc')
        # Get top 10 by p-value
        top_pval_down = down_genes.nlargest(10, '-log10(pvals_adj)')
        # Find intersection of gene names
        down_genes_to_label = set(top_fc_down['names']).intersection(set(top_pval_down['names']))
        # If intersection is empty, take top 5 from each
        if not down_genes_to_label:
            down_genes_to_label = set(top_fc_down['names'].head(5)).union(set(top_pval_down['names'].head(5)))
        
        # Label the selected genes
        for gene_name in down_genes_to_label:
            gene = down_genes[down_genes['names'] == gene_name].iloc[0]
            down_texts.append(plt.text(gene['logfoldchanges'], 
                                gene['-log10(pvals_adj)'],
                                gene['names'],
                                fontsize=18))

    # Adjust text positions to prevent overlap - handle up and down regulated separately first
    if up_texts:
        adjust_text(up_texts,
                arrowprops=dict(arrowstyle='->', color='red', lw=0.5),
                expand_points=(1.8, 1.8),
                force_text=(0.5, 1.0),
                force_points=(0.5, 0.5))
    
    if down_texts:
        adjust_text(down_texts,
                arrowprops=dict(arrowstyle='->', color='blue', lw=0.5),
                expand_points=(1.8, 1.8),
                force_text=(0.5, 1.0),
                force_points=(0.5, 0.5))
    
    # Final adjustment for all texts together to ensure no overlap
    all_texts = up_texts + down_texts
    if all_texts:
        adjust_text(all_texts,
                arrowprops=dict(arrowstyle='->', color='black', lw=0.5),
                expand_points=(2.0, 2.0),
                force_text=(0.7, 1.2),
                force_points=(0.7, 0.7),
                lim=5000)  # Increase iterations for better separation

    # Add axis labels and title
    plt.xlabel('Log2 Fold Change', fontsize=20)
    plt.ylabel('-log10(adjusted p-value)', fontsize=20) 
    plt.title(f'Volcano Plot: {comparator}', fontsize=22)

    # Add horizontal line for p-value cutoff
    plt.axhline(y=-np.log10(0.05), color='k', linestyle='--', alpha=0.3)
    # Add vertical lines for logFC cutoffs
    plt.axvline(x=1, color='k', linestyle='--', alpha=0.3)
    plt.axvline(x=-1, color='k', linestyle='--', alpha=0.3)

    plt.legend(fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    
    os.makedirs('figures/volcano_plots/', exist_ok=True)
    plt.savefig(f'figures/volcano_plots/{comparator}_volcano_plot.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'figures/volcano_plots/{comparator}_volcano_plot.svg', bbox_inches='tight')
    plt.show()

    return deg_df


import gseapy as gp
import time
def pathway_plot(adata, comparator):
    """
    geneset_list = ['KEGG_2021_Human','BioCarta_2016', 'MSigDB_Hallmark_2020', 'Reactome_2022',
                    'CellMarker_2024', 'CellMarker_Augmented_2021']
    """

    results = adata.uns['rank_genes_groups']
    groups = results['names'].dtype.names

    for label in groups:
        gene_set_list = ['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022', 'all']
        color=['darkred', 'darkblue', 'darkgreen']

        # Process both upregulated and downregulated genes
        for direction in ['up', 'down']:
            for num, gene_set in enumerate(gene_set_list):

                if gene_set == 'all':
                    gene_set = ['MSigDB_Hallmark_2020', 'KEGG_2021_Human', 'Reactome_2022']

                groups = results['names'].dtype.names
                dataframes_dict = {}
                for group in groups:
                    df = pd.DataFrame({
                        'names': results['names'][group],
                        'scores': results['scores'][group],
                        'logfoldchanges': results['logfoldchanges'][group],
                        'pvals': results['pvals'][group],
                        'pvals_adj': results['pvals_adj'][group]
                    })
                    dataframes_dict[group] = df

                df_of_interest = dataframes_dict[label]

                print(df_of_interest.head())

                # Invert signs for downregulated analysis
                if direction == 'down':
                    df_of_interest['scores'] = -df_of_interest['scores']
                    df_of_interest['logfoldchanges'] = -df_of_interest['logfoldchanges']

                group_df = df_of_interest
                group_df = group_df[group_df[f'scores'] > 0]
                group_df = group_df[group_df[f'pvals_adj'] < 0.05]
                group_df = group_df[group_df[f'logfoldchanges'] > 0]

                gene_list = group_df['names'].tolist()

                print("Number of genes in the group: ", len(gene_list))

                if len(gene_list) < 30:
                    print(f"No significant pathways found for {gene_set} {direction}regulated")
                    continue

                try:
                    enr = gp.enrichr(gene_list=gene_list, 
                                    gene_sets=gene_set,
                                    organism='Human',
                                    outdir='enrichr_kegg')
                    print(enr.results)

                    # Create figure and plot
                    fig, ax = plt.subplots(figsize=(3, 5))
                    

                    # categorical scatterplot
                    
                    if isinstance(gene_set, list):
                        gp.barplot(enr.results,
                                column="Adjusted P-value",
                                group='Gene_set',
                                size=10,
                                top_term=5,
                                color=color,
                                ax=ax  # Pass the axis object
                                )
                    else:
                        gp.barplot(enr.results,
                                column="Adjusted P-value",
                                group='Gene_set',
                                size=10,
                                top_term=15,
                                color=color[num],
                                ax=ax  # Pass the axis object
                                )
                    
                    # Ensure the output directory exists
                    os.makedirs(f'figures/pathway_{direction}regulated/{comparator}/', exist_ok=True)
                    
                    # Save and close the figure
                    plt.tight_layout()
                    fig.savefig(f"figures/pathway_{direction}regulated/{comparator}/{label}_{gene_set}_barplot.png", 
                               bbox_inches='tight', 
                               dpi=300)
                    fig.savefig(f"figures/pathway_{direction}regulated/{comparator}/{label}_{gene_set}_barplot.svg", 
                               bbox_inches='tight')
                    plt.close(fig)  # Close the figure to free memory
                    
                except ValueError:
                    print(f"No significant pathways found for {comparator} {label} {gene_set} {direction}regulated")
                except Exception as e:
                    print(f"Error for {comparator} {label} {gene_set} {direction}regulated: {e}")
                    # Wait for 5 seconds and then retry
                    time.sleep(5)
                    continue
    return None

def calculate_tissue_proportions(adata, celltypes):
    """Calculate the proportion of each celltype across the entire tissue."""
    # Sum up total proportions for each celltype across all spots
    celltype_totals = {ct: np.sum(adata.obs[ct].values) for ct in celltypes}
    # Calculate total of all celltypes
    total_sum = sum(celltype_totals.values())
    # Calculate relative proportion for each celltype
    return {ct: total/total_sum for ct, total in celltype_totals.items()}


def preprocess_and_validate_data(adata):
    """Preprocess and validate data to avoid numerical issues"""
    logging.info("Preprocessing and validating data...")
    
    # Replace potential infinities and NaNs
    adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)
    
    # Add small epsilon to zeros to avoid log1p issues
    epsilon = 1e-8
    adata.X = adata.X + epsilon
    
    # Basic normalization if not already normalized
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Replace any NaNs that might have been introduced during normalization
    adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)
    
    # Clip extreme values to avoid overflow
    upper_bound = np.percentile(adata.X, 99)  # Use 99th percentile as upper bound
    adata.X = np.clip(adata.X, 0, upper_bound)
    
    # Final NaN check
    adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)
    
    # Run PCA if not already done (needed for Harmony)
    if 'X_pca' not in adata.obsm:
        logging.info("Running PCA...")
        sc.pp.pca(adata)
    
    force_garbage_collection()
    return adata


def batch_correct_adata(combined_adata, batch_key='sample_name'):
    """Apply batch correction to an AnnData object and its layers using Harmony."""
    import scanpy as sc
    import scanpy.external as sce
    import numpy as np
    
    # Store original X matrix
    original_X = combined_adata.X.copy()

    # Add small epsilon to avoid divide by zero
    epsilon = 1e-8
    
    # Check if we have enough samples per batch
    batch_sizes = combined_adata.obs[batch_key].value_counts()
    if (batch_sizes < 2).any():
        print(f"Warning: Some batches have less than 2 samples. Skipping batch correction.")
        return combined_adata

    # First handle the main matrix
    try:
        # Make sure PCA has been run
        if 'X_pca' not in combined_adata.obsm:
            print("Running PCA before Harmony...")
            sc.pp.pca(combined_adata)
        
        # Run Harmony on the PCA matrix
        print(f"Running Harmony on main matrix using {batch_key}...")
        sce.pp.harmony_integrate(combined_adata, batch_key)
        
        # The corrected PCs are now in combined_adata.obsm['X_pca_harmony']
        print("Harmony batch correction completed for main matrix.")
    except Exception as e:
        print(f"Warning: Harmony failed with error {e}. Returning original data.")
        return combined_adata

    # Now handle the layers - we need to correct each layer separately
    for layer_name in combined_adata.layers.keys():
        print(f"Processing layer: {layer_name}")
        try:
            # Create temporary AnnData for this layer
            temp_adata = sc.AnnData(X=combined_adata.layers[layer_name],
                                  var=combined_adata.var,
                                  obs=combined_adata.obs)
            
            # Add small epsilon to avoid divide by zero
            temp_adata.X = temp_adata.X + epsilon
            
            # Normalize and log transform
            sc.pp.normalize_per_cell(temp_adata, counts_per_cell_after=1e4)
            sc.pp.log1p(temp_adata)
            
            # Run PCA and Harmony
            sc.pp.pca(temp_adata)
            sce.pp.harmony_integrate(temp_adata, batch_key)
            
            # Store the harmony-corrected PCs in the layer's obsm
            layer_harmony_key = f'X_pca_harmony_{layer_name}'
            combined_adata.obsm[layer_harmony_key] = temp_adata.obsm['X_pca_harmony']
            
            # Note: We're not modifying the original layer data, just storing the harmony-corrected PCs
            # This is different from the ComBat approach which modified the expression matrix directly
        except Exception as e:
            print(f"Warning: Failed to process layer {layer_name}: {e}")
            continue
    
    return combined_adata

def standardize_layer_name(name):
    """Standardize layer names to be consistent across the codebase."""
    # Keep the original name for the obs column
    return name.replace(' ', '_').replace('(', '').replace(')', '').replace('+', 'plus').replace('-', '_minus')

def standardize_adata(adata, cell_profiles):
    """Standardize an AnnData object's layer names and observation columns to be consistent.
    
    This function:
    1. Creates a mapping between original and standardized cell type names
    2. Standardizes layer names
    3. Adds standardized observation columns while preserving original ones
    4. Returns the modified AnnData object and name mapping
    """
    import logging
    
    # Create mapping between original and standardized names
    name_mapping = {}
    for original_name in cell_profiles.keys():
        std_name = standardize_layer_name(original_name)
        name_mapping[original_name] = std_name
        
    # Standardize layer names
    new_layers = {}
    for key in adata.layers.keys():
        if key.endswith('_genes_pass1'):
            celltype = key.replace('_genes_pass1', '')
            std_key = f"{standardize_layer_name(celltype)}_genes_pass1"
            new_layers[std_key] = adata.layers[key]
            logging.info(f"Standardized layer name from {key} to {std_key}")
    
    # Update layers with standardized names
    adata.layers = new_layers
    
    # Add standardized observation columns while preserving original ones
    for original_name in cell_profiles.keys():
        if original_name in adata.obs:
            std_name = name_mapping[original_name]
            # Only add if the standardized name is different and doesn't already exist
            if std_name != original_name and std_name not in adata.obs:
                adata.obs[std_name] = adata.obs[original_name].copy()
                logging.info(f"Added standardized observation column: {std_name}")
    
    return adata, name_mapping

def standardize_cell_profiles(cell_profiles):
    """Standardize cell profile dictionary keys to match standardized layer names."""
    standardized_profiles = {}
    for key, value in cell_profiles.items():
        std_key = standardize_layer_name(key)
        standardized_profiles[std_key] = value
        logging.info(f"Standardized cell profile key from {key} to {std_key}")
    return standardized_profiles

def process_celltype_data(adata_dict, cell_profiles, output_folder='macrophage_output', license_file=None):
    """Process a single sample's data through the CITEgeist model."""
    from citegeist.core.citegeist_model import CitegeistModel
    import gc
    import logging
    
    processed_adatas = {}
    
    for path, adata in adata_dict.items():
        sample_name = path.split('/')[-2]
        logging.info(f"Processing {sample_name}...")
        
        try:
            # Initialize the model
            model = CitegeistModel(sample_name=sample_name, adata=adata, output_folder=output_folder)
            
            # Load cell profile dictionary
            model.load_cell_profile_dict(cell_profiles)
            
            # Preprocess and run models
            model.split_adata()
            model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=25)
            model.copy_gex_to_protein_adata()
            
            # Preprocess datasets
            model.preprocess_gex()
            model.preprocess_antibody()
            logging.info(f"Model preprocessing complete for {sample_name}")
            
            # Register Gurobi license if provided
            if license_file:
                model.register_gurobi(license_file)
            
            # Plot cell proportions
            model.append_proportions_to_adata(key='finetuned')
            model.append_gex_to_adata(pass_number=1)
            
            prop_gex_adata = model.get_adata()
            
            # Standardize the AnnData object
            prop_gex_adata, _ = standardize_adata(prop_gex_adata, cell_profiles)
            
            processed_adatas[sample_name] = prop_gex_adata
            
            # Clear memory
            del model
            gc.collect()
            
        except Exception as e:
            logging.error(f"Error processing {sample_name}: {e}")
            continue
    
    return processed_adatas

def plot_celltype_proportions(responder_samples, nonresponder_samples, celltypes, chunk_size=5, output_dir='figures/celltype_proportions'):
    """Create comparison plots for celltype proportions between conditions."""
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import stats
    import os
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Split celltypes into chunks
    num_chunks = (len(celltypes) + chunk_size - 1) // chunk_size
    
    for chunk_idx in range(num_chunks):
        chunk_celltypes = celltypes[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size]
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 12))
        x_pos = np.arange(len(chunk_celltypes))
        width = 0.35
        
        # Plot for each group
        for ax, (title, samples) in zip([ax1, ax2], [('Responders', responder_samples), ('Non-Responders', nonresponder_samples)]):
            biopsy_vals = [s['biopsy'] for s in samples]
            surgical_vals = [s['surgical'] for s in samples]
            
            # Calculate means and standard errors
            biopsy_means = [np.mean([b[ct] for b in biopsy_vals]) for ct in chunk_celltypes]
            surgical_means = [np.mean([s[ct] for s in surgical_vals]) for ct in chunk_celltypes]
            biopsy_sems = [stats.sem([b[ct] for b in biopsy_vals]) if len(biopsy_vals) > 1 else 0 for ct in chunk_celltypes]
            surgical_sems = [stats.sem([s[ct] for s in surgical_vals]) if len(surgical_vals) > 1 else 0 for ct in chunk_celltypes]
            
            # Plot bars with error bars
            ax.bar(x_pos - width/2, biopsy_means, width, label='Biopsy', color='skyblue', 
                  yerr=biopsy_sems, capsize=5, error_kw=dict(elinewidth=2, capthick=2))
            ax.bar(x_pos + width/2, surgical_means, width, label='Surgical', color='lightcoral',
                  yerr=surgical_sems, capsize=5, error_kw=dict(elinewidth=2, capthick=2))
            
            # Add significance stars
            for i, ct in enumerate(chunk_celltypes):
                if len(biopsy_vals) > 1 and len(surgical_vals) > 1:
                    _, p_val = stats.ttest_ind([b[ct] for b in biopsy_vals], [s[ct] for s in surgical_vals])
                    height = max(biopsy_means[i], surgical_means[i])
                    max_sem = max(biopsy_sems[i], surgical_sems[i])
                    stars = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''
                    if stars:
                        ax.text(i, height + max_sem + 0.02, stars, ha='center')
            
            ax.set_title(f'{title}: Tissue-wide Celltype Proportions (Biopsy vs Surgical)', fontsize=14)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(chunk_celltypes, rotation=45, ha='right')
            ax.legend()
            ax.grid(True, alpha=0.3)
            ax.set_ylabel('Proportion of Total Tissue')
        
        # Add significance legend
        fig.text(0.02, 0.98, '* p < 0.05\n** p < 0.01\n*** p < 0.001',
                verticalalignment='top', fontsize=14,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/tissue_wide_proportions_chunk_{chunk_idx + 1}.png', 
                    bbox_inches='tight', dpi=300)
        plt.savefig(f'{output_dir}/tissue_wide_proportions_chunk_{chunk_idx + 1}.svg', 
                    bbox_inches='tight')
        plt.close()

def analyze_differential_expression(adata, groupby, output_prefix, condition_key=None):
    """Perform differential expression analysis and create visualizations."""
    import scanpy as sc
    import os
    
    # Create all necessary output directories
    os.makedirs('figures/rank_genes_groups_leiden', exist_ok=True)
    os.makedirs('figures/rank_genes_groups_condition', exist_ok=True)
    os.makedirs('figures/volcano_plots', exist_ok=True)
    os.makedirs('figures/pathway_upregulated', exist_ok=True)
    os.makedirs('figures/pathway_downregulated', exist_ok=True)
    os.makedirs(f'figures/rank_genes_groups_{groupby}', exist_ok=True)
    
    # Perform differential expression analysis
    sc.tl.rank_genes_groups(adata, groupby, method='wilcoxon')
    
    # Generate plots
    sc.pl.rank_genes_groups(adata, save=f"//{output_prefix}_rank_genes_groups.png")
    
    # Create the SVG version using matplotlib's savefig
    # First get the current figure
    fig = plt.gcf()
    # Create directory if not exists (use the same directory pattern as scanpy)
    output_dir = os.path.join('figures', f'rank_genes_groups_{groupby}')
    os.makedirs(output_dir, exist_ok=True)
    # Save as SVG
    fig.savefig(os.path.join(output_dir, f"{output_prefix}_rank_genes_groups.svg"), bbox_inches='tight')
    
    volcano_plot(adata, output_prefix)
    pathway_plot(adata, output_prefix)
    
    return adata

def process_commot_analysis(adata, celltype_profile_dict, distance_threshold=500):
    """Process COMMOT analysis for a single sample."""
    import commot as ct
    
    # Prepare ligand-receptor database
    df_ligrec = ct.pp.ligand_receptor_database(database='CellChat', species='human')
    df_cellchat_filtered = ct.pp.filter_lr_database(df_ligrec, adata, min_cell_pct=0.01)
    
    # Run COMMOT
    ct.tl.spatial_communication(adata, 
                              database_name='cellchat', 
                              df_ligrec=df_cellchat_filtered, 
                              dis_thr=distance_threshold, 
                              heteromeric=True, 
                              pathway_sum=False)
    
    # Extract signals
    sender_signal = adata.obsm['commot-cellchat-sum-sender']
    receiver_signal = adata.obsm['commot-cellchat-sum-receiver']
    
    return {'sender_signal': sender_signal, 'receiver_signal': receiver_signal}


def log_memory_usage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1024 / 1024
    logging.info(f"Memory usage: {mem:.2f} MB")


def calculate_p_value(data1, data2):
    """Calculate p-value between two datasets with explicit handling of missing values."""
    try:
        # If we have raw values, use t-test
        if ('raw_values' in data1 and 'raw_values' in data2 and
            len(data1['raw_values']['group1']) >= 2 and 
            len(data2['raw_values']['group2']) >= 2):
            
            vals1 = data1['raw_values']['group1']
            vals2 = data2['raw_values']['group2']
            
            # Log the values for debugging
            logging.debug(f"Group 1 values: {vals1}")
            logging.debug(f"Group 2 values: {vals2}")
            
            # If all values in one group are 0 and the other has non-zero values,
            # this is likely significant
            if (all(v == 0 for v in vals1) and not all(v == 0 for v in vals2)):
                logging.info("Pathway present only in group 2")
                return 0.01
            elif (all(v == 0 for v in vals2) and not all(v == 0 for v in vals1)):
                logging.info("Pathway present only in group 1")
                return 0.01
            
            # Otherwise perform t-test
                _, p_val = stats.ttest_ind(vals1, vals2)
                return p_val
        
        # If no raw values but we have p-values, use the minimum
        if 'p_val' in data1 and 'p_val' in data2:
            return min(data1['p_val'], data2['p_val'])
        elif 'p_val' in data1:
            return data1['p_val']
        elif 'p_val' in data2:
            return data2['p_val']
        
        # Default to 1.0 (not significant)
        return 1.0
    except Exception as e:
        logging.error(f"Error calculating p-value: {e}")
        return 1.0
    


def create_comparison_plot(celltype, data1, data2, comparison_name, output_prefix):
    """Create comparison plots for a specific cell type between two datasets."""
    try:
        # Add logging at the start
        logging.info(f"Starting comparison for {celltype} in {comparison_name}")
        
        # Get all unique pathways from both groups
        pathways1 = set()
        pathways2 = set()
        
        # Collect all pathways and log counts
        for results in data1.values():
            if celltype in results['sender'].index:
                pathways1.update(col.replace('s-', '') for col in results['sender'].columns)
            if celltype in results['receiver'].index:
                pathways1.update(col.replace('r-', '') for col in results['receiver'].columns)

        for results in data2.values():
            if celltype in results['sender'].index:
                pathways2.update(col.replace('s-', '') for col in results['sender'].columns)
            if celltype in results['receiver'].index:
                pathways2.update(col.replace('r-', '') for col in results['receiver'].columns)
        
        # Use all pathways from either group
        all_pathways = pathways1 | pathways2
        
        # Remove total-total if present
        if 'total-total' in all_pathways:
            all_pathways.remove('total-total')
        
        # Initialize storage for pathway results
        all_pathways_data = {}
        
        # Process each pathway
        for pathway in all_pathways:
            # Collect cell-level values for each group
            group1_sender_cells = []
            group1_receiver_cells = []
            group2_sender_cells = []
            group2_receiver_cells = []
            
            # Collect values from group 1
            for results in data1.values():
                if celltype in results['sender'].index:
                    # Get cell-level values for this pathway
                    cell_mask = results['sender'].index.str.contains(celltype)
                    if f"s-{pathway}" in results['sender'].columns:
                        vals = results['sender'].loc[cell_mask, f"s-{pathway}"].values
                    elif pathway in results['sender'].columns:
                        vals = results['sender'].loc[cell_mask, pathway].values
                    else:
                        vals = np.zeros(sum(cell_mask))
                    group1_sender_cells.extend(vals)
                
                if celltype in results['receiver'].index:
                    cell_mask = results['receiver'].index.str.contains(celltype)
                    if f"r-{pathway}" in results['receiver'].columns:
                        vals = results['receiver'].loc[cell_mask, f"r-{pathway}"].values
                    elif pathway in results['receiver'].columns:
                        vals = results['receiver'].loc[cell_mask, pathway].values
                    else:
                        vals = np.zeros(sum(cell_mask))
                    group1_receiver_cells.extend(vals)
            
            # Collect values from group 2
            for results in data2.values():
                if celltype in results['sender'].index:
                    cell_mask = results['sender'].index.str.contains(celltype)
                    if f"s-{pathway}" in results['sender'].columns:
                        vals = results['sender'].loc[cell_mask, f"s-{pathway}"].values
                    elif pathway in results['sender'].columns:
                        vals = results['sender'].loc[cell_mask, pathway].values
                    else:
                        vals = np.zeros(sum(cell_mask))
                    group2_sender_cells.extend(vals)
                
                if celltype in results['receiver'].index:
                    cell_mask = results['receiver'].index.str.contains(celltype)
                    if f"r-{pathway}" in results['receiver'].columns:
                        vals = results['receiver'].loc[cell_mask, f"r-{pathway}"].values
                    elif pathway in results['receiver'].columns:
                        vals = results['receiver'].loc[cell_mask, pathway].values
                    else:
                        vals = np.zeros(sum(cell_mask))
                    group2_receiver_cells.extend(vals)
            
            # Process sender values
            if group1_sender_cells or group2_sender_cells:
                sender_mean1 = np.mean(group1_sender_cells) if group1_sender_cells else 0.0
                sender_mean2 = np.mean(group2_sender_cells) if group2_sender_cells else 0.0
                
                # Perform Mann-Whitney U test if we have enough samples
                if len(group1_sender_cells) >= 10 and len(group2_sender_cells) >= 10:
                    _, sender_p_val = stats.mannwhitneyu(group1_sender_cells, group2_sender_cells, alternative='two-sided')
                else:
                    sender_p_val = 1.0
                
                all_pathways_data[pathway] = all_pathways_data.get(pathway, {})
                all_pathways_data[pathway]['sender'] = {
                    'group1_mean': sender_mean1,
                    'group2_mean': sender_mean2,
                    'p_val': sender_p_val,
                    'n_cells_1': len(group1_sender_cells),
                    'n_cells_2': len(group2_sender_cells)
                }
            
            # Process receiver values
            if group1_receiver_cells or group2_receiver_cells:
                receiver_mean1 = np.mean(group1_receiver_cells) if group1_receiver_cells else 0.0
                receiver_mean2 = np.mean(group2_receiver_cells) if group2_receiver_cells else 0.0
                
                # Perform Mann-Whitney U test if we have enough samples
                if len(group1_receiver_cells) >= 10 and len(group2_receiver_cells) >= 10:
                    _, receiver_p_val = stats.mannwhitneyu(group1_receiver_cells, group2_receiver_cells, alternative='two-sided')
                else:
                    receiver_p_val = 1.0
                
                all_pathways_data[pathway] = all_pathways_data.get(pathway, {})
                all_pathways_data[pathway]['receiver'] = {
                    'group1_mean': receiver_mean1,
                    'group2_mean': receiver_mean2,
                    'p_val': receiver_p_val,
                    'n_cells_1': len(group1_receiver_cells),
                    'n_cells_2': len(group2_receiver_cells)
                }
        
        # Log significant findings
        sig_count = sum(1 for p_data in all_pathways_data.values()
                       if min(p_data.get('sender', {}).get('p_val', 1.0),
                            p_data.get('receiver', {}).get('p_val', 1.0)) < 0.05)
        logging.info(f"Found {sig_count} significant pathways for {celltype}")
        
        return all_pathways_data
        
    except Exception as e:
        logging.error(f"Error in create_comparison_plot for {celltype}: {e}")
        logging.exception("Detailed traceback:")
        return None
    



def plot_dot_plot(ax, sender_data, receiver_data, pathway_name, group_name):
    """Create a dot plot for a specific pathway showing sender-receiver relationships.
    
    Args:
        ax: matplotlib axis
        sender_data: DataFrame with sender signals for the pathway
        receiver_data: DataFrame with receiver signals for the pathway
        pathway_name: name of the pathway being plotted (without s- or r- prefix)
        group_name: name of the group for the plot title
    """
    # Input validation
    assert ax is not None, "Matplotlib axis is None"
    assert isinstance(sender_data, pd.DataFrame), f"sender_data must be a DataFrame, got {type(sender_data)}"
    assert isinstance(receiver_data, pd.DataFrame), f"receiver_data must be a DataFrame, got {type(receiver_data)}"
    
    # Add prefixes for accessing the correct columns
    sender_col = pathway_name if pathway_name.startswith('s-') else f"s-{pathway_name}"
    receiver_col = pathway_name if pathway_name.startswith('r-') else f"r-{pathway_name}"
    
    logging.info(f"Starting dot plot for pathway {pathway_name} in {group_name}")
    logging.info(f"Looking for sender column {sender_col} and receiver column {receiver_col}")
    
    # Verify columns exist
    if sender_col not in sender_data.columns or receiver_col not in receiver_data.columns:
        # Try without prefixes if not found
        if pathway_name in sender_data.columns and pathway_name in receiver_data.columns:
            sender_col = pathway_name
            receiver_col = pathway_name
        else:
            logging.warning(f"Missing columns for pathway {pathway_name}. Sender: {sender_col in sender_data.columns}, Receiver: {receiver_col in receiver_data.columns}")
            return False
    
    # Get cell types (should be the same for both sender and receiver)
    celltypes = sorted(list(set(sender_data.index) | set(receiver_data.index)))
    assert len(celltypes) > 0, "No cell types found in data"
    logging.info(f"Cell types found: {celltypes}")
    
    # Create matrix for sender-receiver interactions
    n_cells = len(celltypes)
    signal_matrix = np.zeros((n_cells, n_cells))
    
    # Fill the matrix
    for i, sender in enumerate(celltypes):
        for j, receiver in enumerate(celltypes):
            if sender in sender_data.index and receiver in receiver_data.index:
                sender_val = sender_data.loc[sender, sender_col]
                receiver_val = receiver_data.loc[receiver, receiver_col]
                if sender_val > 0 and receiver_val > 0:
                    signal_matrix[i, j] = np.sqrt(sender_val * receiver_val)
                    logging.info(f"Matrix position [{i}, {j}] filled with value {signal_matrix[i, j]} for {sender}->{receiver}")
    
    # Verify we have some non-zero values
    if not np.any(signal_matrix > 0):
        logging.warning(f"No valid signals found for pathway {pathway_name}")
        return False
    
    # Create the dot plot
    max_value = np.max(signal_matrix)
    scatter_size = 1000  # Base size
    logging.info(f"Maximum signal value: {max_value}")
    
    # Plot dots with raw values (remove normalization)
    dots_plotted = 0
    for i in range(n_cells):
        for j in range(n_cells):
            if signal_matrix[i, j] > 0:
                size = scatter_size * signal_matrix[i, j]  # Use raw value for size
                color = plt.cm.viridis(0.3 + 0.7 * signal_matrix[i, j])  # Use raw value for color
                ax.scatter(j, i, s=size, color=color, edgecolor='black', linewidth=1, alpha=0.8)
                dots_plotted += 1
    
    if dots_plotted == 0:
        logging.warning(f"No dots plotted for pathway {pathway_name}")
        return False
        
    logging.info(f"Total dots plotted: {dots_plotted}")
    
    # Customize plot
    ax.set_xticks(range(n_cells))
    ax.set_yticks(range(n_cells))
    ax.set_xticklabels(celltypes, rotation=45, ha='right', fontsize=14)
    ax.set_yticklabels(celltypes, fontsize=14)
    ax.set_xlabel('Receiver Cell Types', fontsize=14)
    ax.set_ylabel('Sender Cell Types', fontsize=14)
    ax.set_title(f'{group_name}\n{pathway_name}', fontsize=16)
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Update colorbar to show raw values
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis)
    sm.set_array([0, np.max(signal_matrix)])  # Use actual range of values
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.05, aspect=40)
    cbar.set_label('Raw Signal Strength', fontsize=14)
    
    return True

def find_significant_pathways_by_comparison(commot_results, comparison):
    """Find significant pathways for a specific comparison on a celltype by celltype basis.
    
    Args:
        commot_results (dict): Dictionary of raw COMMOT results by sample
        comparison (tuple): Tuple containing (group1_key, group2_key, comparison_name, output_prefix)
    
    Returns:
        dict: Dictionary of significant pathways by celltype
    """
    def extract_celltype(index_name):
        # Remove '_genes_pass1' suffix if present
        if '_genes_pass1' in index_name:
            index_name = index_name.replace('_genes_pass1', '')
        # Remove barcode prefix (everything before and including first '_')
        if '_' in index_name:
            index_name = index_name.split('_', 1)[1]
        return index_name
    
    # Extract comparison information
    group1_key, group2_key, comparison_name, output_prefix = comparison
    
    # Separate samples into groups
    group1_samples = commot_results[group1_key]
    group2_samples = commot_results[group2_key]
    
    print(group1_samples)
    print(group2_samples)
    # Initialize results dictionary
    significant_results = {
        'sender': {},
        'receiver': {}
    }
    
    # Get all unique pathways from both groups
    all_pathways = set()
    for data_dict in [group1_samples, group2_samples]:
        for sample_data in data_dict.values():
            sender_paths = {col.replace('s-', '') for col in sample_data['sender_signal'].columns}
            receiver_paths = {col.replace('r-', '') for col in sample_data['receiver_signal'].columns}
            all_pathways.update(sender_paths | receiver_paths)
    
    # Remove total-total if present
    if 'total-total' in all_pathways:
        all_pathways.remove('total-total')
    
    # Get all unique celltypes
    all_celltypes = set()
    for data_dict in [group1_samples, group2_samples]:
        for sample_data in data_dict.values():
            sender_types = set(sample_data['sender_signal'].index.map(extract_celltype))
            receiver_types = set(sample_data['receiver_signal'].index.map(extract_celltype))
            all_celltypes.update(sender_types | receiver_types)
    
    # For each celltype
    for celltype in all_celltypes:
        significant_results['sender'][celltype] = {}
        significant_results['receiver'][celltype] = {}
        
        # For each pathway
        for pathway in all_pathways:
            # Initialize lists for values
            group1_sender_vals = []
            group1_receiver_vals = []
            group2_sender_vals = []
            group2_receiver_vals = []
            
            # Collect values from group 1
            for sample_data in group1_samples.values():
                sender_col = f"s-{pathway}"
                receiver_col = f"r-{pathway}"
                
                # Get sender values for spots of this celltype
                sender_mask = sample_data['sender_signal'].index.map(extract_celltype) == celltype
                if sender_col in sample_data['sender_signal'].columns:
                    vals = sample_data['sender_signal'].loc[sender_mask, sender_col].values
                    group1_sender_vals.extend(vals)
                
                # Get receiver values for spots of this celltype
                receiver_mask = sample_data['receiver_signal'].index.map(extract_celltype) == celltype
                if receiver_col in sample_data['receiver_signal'].columns:
                    vals = sample_data['receiver_signal'].loc[receiver_mask, receiver_col].values
                    group1_receiver_vals.extend(vals)
            
            # Collect values from group 2
            for sample_data in group2_samples.values():
                sender_col = f"s-{pathway}"
                receiver_col = f"r-{pathway}"
                
                # Get sender values for spots of this celltype
                sender_mask = sample_data['sender_signal'].index.map(extract_celltype) == celltype
                if sender_col in sample_data['sender_signal'].columns:
                    vals = sample_data['sender_signal'].loc[sender_mask, sender_col].values
                    group2_sender_vals.extend(vals)
                
                # Get receiver values for spots of this celltype
                receiver_mask = sample_data['receiver_signal'].index.map(extract_celltype) == celltype
                if receiver_col in sample_data['receiver_signal'].columns:
                    vals = sample_data['receiver_signal'].loc[receiver_mask, receiver_col].values
                    group2_receiver_vals.extend(vals)
            
            # Add zeros for missing pathways
            if not group1_sender_vals and group2_sender_vals:
                group1_sender_vals = [0.0] * len(group2_sender_vals)
            if not group2_sender_vals and group1_sender_vals:
                group2_sender_vals = [0.0] * len(group1_sender_vals)
            if not group1_receiver_vals and group2_receiver_vals:
                group1_receiver_vals = [0.0] * len(group2_receiver_vals)
            if not group2_receiver_vals and group1_receiver_vals:
                group2_receiver_vals = [0.0] * len(group1_receiver_vals)
            
            # Perform statistical tests if we have enough samples
            if len(group1_sender_vals) >= 1 and len(group2_sender_vals) >= 1:
                try:
                    # Convert values to numpy arrays
                    group1_sender_arr = np.array(group1_sender_vals, dtype=float)
                    group2_sender_arr = np.array(group2_sender_vals, dtype=float)
                    
                    # T-test for sender values
                    _, p_val = stats.ttest_ind(group1_sender_arr, group2_sender_arr)
                    if not np.isnan(p_val):
                        significant_results['sender'][celltype][pathway] = {
                            'p_value': p_val,
                            'group1_mean': np.mean(group1_sender_arr),
                            'group2_mean': np.mean(group2_sender_arr),
                            'log2fc': np.log2(np.mean(group2_sender_arr + 1) / np.mean(group1_sender_arr + 1))
                        }
                except Exception as e:
                    logging.debug(f"Error in sender t-test for {pathway}, {celltype}: {e}")
            
            if len(group1_receiver_vals) >= 1 and len(group2_receiver_vals) >= 1:
                try:
                    # Convert values to numpy arrays
                    group1_receiver_arr = np.array(group1_receiver_vals, dtype=float) 
                    group2_receiver_arr = np.array(group2_receiver_vals, dtype=float)
                    
                    # T-test for receiver values
                    _, p_val = stats.ttest_ind(group1_receiver_arr, group2_receiver_arr)
                    if not np.isnan(p_val):
                        significant_results['receiver'][celltype][pathway] = {
                            'p_value': p_val,
                            'group1_mean': np.mean(group1_receiver_arr),
                            'group2_mean': np.mean(group2_receiver_arr),
                            'log2fc': np.log2(np.mean(group2_receiver_arr + 1) / np.mean(group1_receiver_arr + 1))
                        }
                except Exception as e:
                    logging.debug(f"Error in receiver t-test for {pathway}, {celltype}: {e}")

    # Apply multiple test correction after collecting all p-values
    from statsmodels.stats.multitest import multipletests
    
    # Collect all p-values and their corresponding results
    all_pvals = []
    all_results = []
    
    # Collect sender p-values
    for celltype in significant_results['sender']:
        for pathway in significant_results['sender'][celltype]:
            result = significant_results['sender'][celltype][pathway]
            all_pvals.append(result['p_value'])
            all_results.append(('sender', celltype, pathway, result))
    
    # Collect receiver p-values
    for celltype in significant_results['receiver']:
        for pathway in significant_results['receiver'][celltype]:
            result = significant_results['receiver'][celltype][pathway]
            all_pvals.append(result['p_value'])
            all_results.append(('receiver', celltype, pathway, result))
    
    if all_pvals:
        # Apply Benjamini-Hochberg correction
        _, pvals_corrected, _, _ = multipletests(all_pvals, method='fdr_bh')
        
        # Update results with corrected p-values
        for (signal_type, celltype, pathway, result), corrected_pval in zip(all_results, pvals_corrected):
            result['p_value_corrected'] = corrected_pval
            # Only keep significant results after correction
            if corrected_pval >= 0.05:
                if signal_type == 'sender':
                    del significant_results['sender'][celltype][pathway]
                else:
                    del significant_results['receiver'][celltype][pathway]
    
    logging.info(f"Significant results: {significant_results}")

    return significant_results

def plot_significant_results_by_celltype(significant_results, comparison, output_dir='figures/signaling/celltype_comparisons'):
    """
    Create bar plots and dot plots for each celltype showing significant sender and receiver signals.
    
    Args:
        significant_results (dict): Dictionary containing significant results
        comparison (tuple): Tuple containing (group1_key, group2_key, comparison_name, output_prefix)
        output_dir (str): Base directory to save the plots
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    
    # Unpack comparison tuple
    group1_key, group2_key, comparison_name, output_prefix = comparison
    
    # Create specific output directory for this comparison
    comparison_dir = os.path.join(output_dir, output_prefix)
    os.makedirs(comparison_dir, exist_ok=True)
    
    # Get all unique celltypes
    all_celltypes = set(significant_results['sender'].keys()) | set(significant_results['receiver'].keys())
    
    for celltype in all_celltypes:
        # Get significant pathways for this celltype
        sender_pathways = significant_results['sender'].get(celltype, {})
        receiver_pathways = significant_results['receiver'].get(celltype, {})
        
        if not sender_pathways and not receiver_pathways:
            continue
            
        # Combine all pathways with significant results
        all_significant_pathways = set(sender_pathways.keys()) | set(receiver_pathways.keys())
        
        if not all_significant_pathways:
            continue
            
        # Create Bar Plot figure
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Set up positions for bars
        width = 0.35
        
        # --- Plot Sender Signals (Bar Plot) ---
        sender_group1_means = []
        sender_group2_means = []
        sender_pathway_labels = []
        
        for pathway in sorted(all_significant_pathways):  # Sort for consistency
            if pathway in sender_pathways:
                pathway_data = sender_pathways[pathway]
                sender_group1_means.append(pathway_data['group1_mean'])
                sender_group2_means.append(pathway_data['group2_mean'])
                sender_pathway_labels.append(pathway)
        
        if sender_pathway_labels: # Check if there's data to plot
            x_sender = np.arange(len(sender_pathway_labels))
            
            ax1.bar(x_sender - width/2, sender_group1_means, width, 
                   label=group1_key.replace('_', ' ').title(), 
                   color='skyblue')
            ax1.bar(x_sender + width/2, sender_group2_means, width, 
                   label=group2_key.replace('_', ' ').title(), 
                   color='lightcoral')
            
            # Calculate the maximum height for setting y-axis limits
            max_height = max(max(sender_group1_means or [0]), max(sender_group2_means or [0]))
            
            # Add significance stars using x_sender
            for i, pathway in enumerate(sender_pathway_labels):
                p_val = sender_pathways[pathway]['p_value_corrected']
                height = max(sender_group1_means[i], sender_group2_means[i])
                stars = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''
                if stars:
                    ax1.text(x_sender[i], height * 1.05, stars, ha='center')
            
            # Set y-axis limit with extra space for stars
            ax1.set_ylim(0, max_height * 1.25)  # 25% extra space above the highest bar
            
            ax1.set_xticks(x_sender)
            ax1.set_xticklabels(sender_pathway_labels, rotation=45, ha='right')
            ax1.set_title(f'{celltype} - Sender Signals (Bar Plot)\n{comparison_name}')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            ax1.set_ylabel("Mean Signal")

            # --- Create Sender Dot Plot ---
            fig_sender_dot, ax_sender_dot = plt.subplots(figsize=(8, max(6, len(sender_pathway_labels) * 0.3)))
            x_dot_labels = [group1_key.replace('_', ' ').title(), group2_key.replace('_', ' ').title()]
            y_dot_labels = sender_pathway_labels

            # Calculate dot sizes
            all_sender_means = sender_group1_means + sender_group2_means
            max_sender_mean = max(all_sender_means) if all_sender_means else 0
            base_size = 20  # Minimum dot size
            scale_factor = 500 / max_sender_mean if max_sender_mean > 0 else 0
            sizes1 = [base_size + m * scale_factor for m in sender_group1_means]
            sizes2 = [base_size + m * scale_factor for m in sender_group2_means]
            
            y_pos_dot = np.arange(len(y_dot_labels))
            ax_sender_dot.scatter([0] * len(y_dot_labels), y_pos_dot, s=sizes1, color='skyblue', alpha=0.7, label=x_dot_labels[0])
            ax_sender_dot.scatter([1] * len(y_dot_labels), y_pos_dot, s=sizes2, color='lightcoral', alpha=0.7, label=x_dot_labels[1])

            # Add size legend (example sizes)
            legend_sizes = [max_sender_mean * 0.1, max_sender_mean * 0.5, max_sender_mean] if max_sender_mean > 0 else [0.1, 0.5, 1.0]
            legend_markers = [plt.scatter([], [], s=base_size + s * scale_factor, color='gray', alpha=0.7) for s in legend_sizes]
            ax_sender_dot.legend(legend_markers, [f'{s:.2f}' for s in legend_sizes], title='Mean Signal', loc='upper right', bbox_to_anchor=(1.25, 1))


            ax_sender_dot.set_xticks([0, 1])
            ax_sender_dot.set_xticklabels(x_dot_labels, rotation=45, ha='right')
            ax_sender_dot.set_yticks(y_pos_dot)
            ax_sender_dot.set_yticklabels(y_dot_labels)
            ax_sender_dot.set_title(f'{celltype} - Sender Signals (Dot Plot)\n{comparison_name}')
            ax_sender_dot.invert_yaxis()
            ax_sender_dot.grid(True, axis='x', linestyle='--', alpha=0.5)

            safe_celltype = celltype.replace('/', '_').replace(' ', '_')
            dot_filename = f'{comparison_dir}/{safe_celltype}_sender_dotplot.png'
            fig_sender_dot.savefig(dot_filename, bbox_inches='tight', dpi=300)
            fig_sender_dot.savefig(dot_filename.replace('.png', '.svg'), bbox_inches='tight')
            plt.close(fig_sender_dot)

        else:
            # If no sender pathways, add placeholder text to bar plot axis
             ax1.text(0.5, 0.5, 'No significant sender pathways', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
             ax1.set_title(f'{celltype} - Sender Signals (Bar Plot)\n{comparison_name}')


        # --- Plot Receiver Signals (Bar Plot) ---
        receiver_group1_means = []
        receiver_group2_means = []
        receiver_pathway_labels = []
        
        for pathway in sorted(all_significant_pathways):  # Sort for consistency
            if pathway in receiver_pathways:
                pathway_data = receiver_pathways[pathway]
                receiver_group1_means.append(pathway_data['group1_mean'])
                receiver_group2_means.append(pathway_data['group2_mean'])
                receiver_pathway_labels.append(pathway)
        
        if receiver_pathway_labels: # Check if there's data to plot
            x_receiver = np.arange(len(receiver_pathway_labels))
            
            ax2.bar(x_receiver - width/2, receiver_group1_means, width,
                   label=group1_key.replace('_', ' ').title(),
                   color='skyblue')
            ax2.bar(x_receiver + width/2, receiver_group2_means, width,
                   label=group2_key.replace('_', ' ').title(),
                   color='lightcoral')
            
            # Calculate the maximum height for setting y-axis limits
            max_height = max(max(receiver_group1_means or [0]), max(receiver_group2_means or [0]))
            
            # Add significance stars using x_receiver
            for i, pathway in enumerate(receiver_pathway_labels):
                p_val = receiver_pathways[pathway]['p_value_corrected']
                height = max(receiver_group1_means[i], receiver_group2_means[i])
                stars = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''
                if stars:
                    ax2.text(x_receiver[i], height * 1.05, stars, ha='center')
            
            # Set y-axis limit with extra space for stars
            ax2.set_ylim(0, max_height * 1.25)  # 25% extra space above the highest bar

            ax2.set_xticks(x_receiver)
            ax2.set_xticklabels(receiver_pathway_labels, rotation=45, ha='right')
            ax2.set_title(f'{celltype} - Receiver Signals (Bar Plot)\n{comparison_name}')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            ax2.set_ylabel("Mean Signal")


            # --- Create Receiver Dot Plot ---
            fig_receiver_dot, ax_receiver_dot = plt.subplots(figsize=(8, max(6, len(receiver_pathway_labels) * 0.3)))
            x_dot_labels = [group1_key.replace('_', ' ').title(), group2_key.replace('_', ' ').title()]
            y_dot_labels = receiver_pathway_labels

            # Calculate dot sizes
            all_receiver_means = receiver_group1_means + receiver_group2_means
            max_receiver_mean = max(all_receiver_means) if all_receiver_means else 0
            # base_size = 20 (already defined)
            scale_factor = 500 / max_receiver_mean if max_receiver_mean > 0 else 0
            sizes1 = [base_size + m * scale_factor for m in receiver_group1_means]
            sizes2 = [base_size + m * scale_factor for m in receiver_group2_means]

            y_pos_dot = np.arange(len(y_dot_labels))
            ax_receiver_dot.scatter([0] * len(y_dot_labels), y_pos_dot, s=sizes1, color='skyblue', alpha=0.7, label=x_dot_labels[0])
            ax_receiver_dot.scatter([1] * len(y_dot_labels), y_pos_dot, s=sizes2, color='lightcoral', alpha=0.7, label=x_dot_labels[1])

            # Add size legend (example sizes)
            legend_sizes = [max_receiver_mean * 0.1, max_receiver_mean * 0.5, max_receiver_mean] if max_receiver_mean > 0 else [0.1, 0.5, 1.0]
            legend_markers = [plt.scatter([], [], s=base_size + s * scale_factor, color='gray', alpha=0.7) for s in legend_sizes]
            ax_receiver_dot.legend(legend_markers, [f'{s:.2f}' for s in legend_sizes], title='Mean Signal', loc='upper right', bbox_to_anchor=(1.25, 1))

            ax_receiver_dot.set_xticks([0, 1])
            ax_receiver_dot.set_xticklabels(x_dot_labels, rotation=45, ha='right')
            ax_receiver_dot.set_yticks(y_pos_dot)
            ax_receiver_dot.set_yticklabels(y_dot_labels)
            ax_receiver_dot.set_title(f'{celltype} - Receiver Signals (Dot Plot)\n{comparison_name}')
            ax_receiver_dot.invert_yaxis()
            ax_receiver_dot.grid(True, axis='x', linestyle='--', alpha=0.5)

            safe_celltype = celltype.replace('/', '_').replace(' ', '_')
            dot_filename = f'{comparison_dir}/{safe_celltype}_receiver_dotplot.png'
            fig_receiver_dot.savefig(dot_filename, bbox_inches='tight', dpi=300)
            fig_receiver_dot.savefig(dot_filename.replace('.png', '.svg'), bbox_inches='tight')
            plt.close(fig_receiver_dot)

        else:
            # If no receiver pathways, add placeholder text to bar plot axis
            ax2.text(0.5, 0.5, 'No significant receiver pathways', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
            ax2.set_title(f'{celltype} - Receiver Signals (Bar Plot)\n{comparison_name}')
        
        # --- Finalize Bar Plot ---
        # Add significance legend to the original bar plot figure
        fig.text(0.02, 0.98, f'* p < 0.05\\n** p < 0.01\\n*** p < 0.001',
                verticalalignment='top', fontsize=14,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        fig.tight_layout(rect=[0, 0, 1, 0.98]) # Adjust layout to prevent overlap with title/legend
        
        # Save the bar plot figure
        safe_celltype = celltype.replace('/', '_').replace(' ', '_')
        bar_filename = f'{comparison_dir}/{safe_celltype}_bar_signals.png' # Rename to avoid overwriting if script is rerun
        fig.savefig(bar_filename, bbox_inches='tight', dpi=300)
        fig.savefig(bar_filename.replace('.png', '.svg'), bbox_inches='tight')
        plt.close(fig) # Close the bar plot figure


def find_significant_pathways(sample_categories, all_celltypes):
    """Find significant pathways for each comparison."""
    significant_pathways_by_comparison = {}
    
    # Print sample structure for debugging
    logging.debug("\nDEBUG: Sample Categories Structure:")
    for category, samples in sample_categories.items():
        logging.debug(f"\n{category}:")
        logging.debug(f"Number of samples: {len(samples)}")
        if samples:  # If we have any samples
            # Print structure of first sample
            first_sample = next(iter(samples.values()))
            logging.debug("\nExample sample structure:")
            logging.debug("Sender shape:", first_sample['sender'].shape)
            logging.debug("Receiver shape:", first_sample['receiver'].shape)
            logging.debug("\nFirst few sender columns:", list(first_sample['sender'].columns)[:5])
            logging.debug("First few receiver columns:", list(first_sample['receiver'].columns)[:5])
            logging.debug("Sender index (cell types):", list(first_sample['sender'].index))
            logging.debug("Receiver index (cell types):", list(first_sample['receiver'].index))
    
    logging.debug("\nDEBUG: All celltypes being used:", sorted(all_celltypes))
    
    # Define the comparisons
    comparisons = [
        ('non_progressing_biopsy', 'non_progressing_surgical'),
        ('progressing_biopsy', 'progressing_surgical'),
        ('non_progressing_biopsy', 'progressing_biopsy'),
        ('non_progressing_surgical', 'progressing_surgical')
    ]
    
    # For each comparison
    for group1_key, group2_key in comparisons:
        comparison_key = f"{group1_key}_vs_{group2_key}"
        significant_pathways_by_comparison[comparison_key] = set()
        
        if group1_key in sample_categories and group2_key in sample_categories:
            group1_data = sample_categories[group1_key]
            group2_data = sample_categories[group2_key]
            
            # Get all pathways from both groups
            pathways1 = set()
            pathways2 = set()
            
            # Extract pathways from first group
            for sample_data in group1_data.values():
                for col in sample_data['sender'].columns:
                    pathways1.add(col[2:] if col.startswith('s-') else col)
                for col in sample_data['receiver'].columns:
                    pathways1.add(col[2:] if col.startswith('r-') else col)
            
            # Extract pathways from second group
            for sample_data in group2_data.values():
                for col in sample_data['sender'].columns:
                    pathways2.add(col[2:] if col.startswith('s-') else col)
                for col in sample_data['receiver'].columns:
                    pathways2.add(col[2:] if col.startswith('r-') else col)
            
            # Use all pathways from either group
            all_pathways = pathways1 | pathways2
            if 'total-total' in all_pathways:
                all_pathways.remove('total-total')
            
            logging.debug(f"\nDEBUG: Comparison {group1_key} vs {group2_key}:")
            logging.debug(f"  Group 1 samples: {len(group1_data)}")
            logging.debug(f"  Group 2 samples: {len(group2_data)}")
            logging.debug(f"  Group 1 pathways: {len(pathways1)}")
            logging.debug(f"  Group 2 pathways: {len(pathways2)}")
            logging.debug(f"  Total unique pathways: {len(all_pathways)}")
            logging.debug(f"  Example pathways: {list(all_pathways)[:5]}")
            
            # For each pathway and celltype combination
            for pathway in all_pathways:
                for celltype in all_celltypes:
                    # Initialize value lists
                    group1_sender_vals = []
                    group1_receiver_vals = []
                    group2_sender_vals = []
                    group2_receiver_vals = []
                    
                    # Collect values from group 1
                    for sample_data in group1_data.values():
                        if celltype in sample_data['sender'].index:
                            sender_col = f"s-{pathway}" if f"s-{pathway}" in sample_data['sender'].columns else pathway
                            if sender_col in sample_data['sender'].columns:
                                val = sample_data['sender'].loc[celltype, sender_col]
                                group1_sender_vals.append(val)
                            else:
                                group1_sender_vals.append(0.0)  # Missing pathway counts as zero
                                
                        if celltype in sample_data['receiver'].index:
                            receiver_col = f"r-{pathway}" if f"r-{pathway}" in sample_data['receiver'].columns else pathway
                            if receiver_col in sample_data['receiver'].columns:
                                val = sample_data['receiver'].loc[celltype, receiver_col]
                                group1_receiver_vals.append(val)
                            else:
                                group1_receiver_vals.append(0.0)  # Missing pathway counts as zero
                    
                    # Collect values from group 2 (same logic)
                    for sample_data in group2_data.values():
                        if celltype in sample_data['sender'].index:
                            sender_col = f"s-{pathway}" if f"s-{pathway}" in sample_data['sender'].columns else pathway
                            if sender_col in sample_data['sender'].columns:
                                val = sample_data['sender'].loc[celltype, sender_col]
                                group2_sender_vals.append(val)
                            else:
                                group2_sender_vals.append(0.0)  # Missing pathway counts as zero
                                
                        if celltype in sample_data['receiver'].index:
                            receiver_col = f"r-{pathway}" if f"r-{pathway}" in sample_data['receiver'].columns else pathway
                            if receiver_col in sample_data['receiver'].columns:
                                val = sample_data['receiver'].loc[celltype, receiver_col]
                                group2_receiver_vals.append(val)
                            else:
                                group2_receiver_vals.append(0.0)  # Missing pathway counts as zero
                    
                    # Print some debugging info for the first few pathways
                    if pathway in list(all_pathways)[:2]:
                        logging.debug(f"\nDEBUG: Values for pathway {pathway}, celltype {celltype}")
                        logging.debug(f"Group 1 sender values: {group1_sender_vals}")
                        logging.debug(f"Group 2 sender values: {group2_sender_vals}")
                        logging.debug(f"Group 1 receiver values: {group1_receiver_vals}")
                        logging.debug(f"Group 2 receiver values: {group2_receiver_vals}")
                    
                    
                    # Only test if we have enough samples
                    if len(group1_sender_vals) >= 1 and len(group2_sender_vals) >= 1:
                        # Perform t-test for sender values
                        try:
                            _, p_val = stats.ttest_ind(group1_sender_vals, group2_sender_vals)
                            if p_val < 0.05 and not np.isnan(p_val):
                                significant_pathways_by_comparison[comparison_key].add(pathway)
                                logging.debug(f"\nDEBUG: Found significant pathway (sender t-test): {pathway} in {celltype}")
                                logging.debug(f"p-value: {p_val:.4f}")
                                logging.debug(f"Group 1 values: {group1_sender_vals}")
                                logging.debug(f"Group 2 values: {group2_sender_vals}")
                        except Exception as e:
                            logging.debug(f"Error in sender t-test for {pathway}, {celltype}: {e}")
                    
                    # Same for receiver values
                    if len(group1_receiver_vals) >= 2 and len(group2_receiver_vals) >= 2:
                        try:
                            _, p_val = stats.ttest_ind(group1_receiver_vals, group2_receiver_vals)
                            if p_val < 0.05 and not np.isnan(p_val):
                                significant_pathways_by_comparison[comparison_key].add(pathway)
                                logging.debug(f"\nDEBUG: Found significant pathway (receiver t-test): {pathway} in {celltype}")
                                logging.debug(f"p-value: {p_val:.4f}")
                                logging.debug(f"Group 1 values: {group1_receiver_vals}")
                                logging.debug(f"Group 2 values: {group2_receiver_vals}")
                        except Exception as e:
                            logging.debug(f"Error in receiver t-test for {pathway}, {celltype}: {e}")
    
    return significant_pathways_by_comparison

def create_pathway_dot_plots(sample_categories, significant_pathways):
    """Create dot plots for each significant pathway for each group."""
    try:
        # Close any existing figures at the start
        plt.close('all')
        
        # For each group
        for group_name, group_data in sample_categories.items():
            logging.info(f"Processing group: {group_name}")
            
            # Average the data across samples in the group
            avg_sender = pd.concat([d['sender'] for d in group_data.values()]).groupby(level=0).mean()
            avg_receiver = pd.concat([d['receiver'] for d in group_data.values()]).groupby(level=0).mean()

            print(significant_pathways)
            
            # Create a dot plot for each significant pathway
            for pathway in sorted(significant_pathways):
                print(pathway)
                print(avg_sender.columns)
                print(avg_receiver.columns)
                if f"s-{pathway}" in avg_sender.columns and f"r-{pathway}" in avg_receiver.columns:
                    logging.info(f"Creating plot for pathway: {pathway}")
                    
                    # Create figure
                    fig, ax = plt.subplots(figsize=(12, 10))
                    
                    # Create the plot
                    success = plot_dot_plot(ax, avg_sender, avg_receiver, pathway, group_name)
                    
                    if success:
                        # Save the plot
                        plt.tight_layout()
                        safe_pathway = pathway.replace("/", "_").replace(" ", "_")
                        safe_group = group_name.replace(" ", "_")
                        filename = f'figures/signaling/{safe_pathway}_{safe_group}_dotplot.png'
                        fig.savefig(filename, bbox_inches='tight', dpi=300)
                        fig.savefig(filename.replace('.png', '.svg'), bbox_inches='tight')
                        logging.info(f"Saved dot plot: {filename}")
                    
                    plt.close(fig)
        
    except Exception as e:
        logging.error(f"Error in create_pathway_dot_plots: {e}")
        logging.exception("Detailed traceback:")
        plt.close('all')

def create_split_signal_plots(sample_categories, significant_pathways):
    """Create split signal plots for each group using only significant pathways for that comparison."""
    try:
        # Close any existing figures at the start
        plt.close('all')
        
        # Initialize storage for all values to find global maximums
        all_sender_values = []
        all_receiver_values = []
        
        # Dictionary to map pathways to their column names
        pathway_to_sender_col = {}
        pathway_to_receiver_col = {}
        
        # First pass to collect all values and map pathways to columns
        for group_data in sample_categories.values():
            for sample_data in group_data.values():
                # Add sender values
                sender_df = sample_data['sender']
                all_sender_values.append(sender_df)
                
                # Map sender pathways to columns
                for col in sender_df.columns:
                    if col.startswith('s-'):
                        pathway = col[2:]  # Remove s- prefix
                        pathway_to_sender_col[pathway] = col
                    else:
                        pathway_to_sender_col[col] = col
                
                # Add receiver values
                receiver_df = sample_data['receiver']
                all_receiver_values.append(receiver_df)
                
                # Map receiver pathways to columns
                for col in receiver_df.columns:
                    if col.startswith('r-'):
                        pathway = col[2:]  # Remove r- prefix
                        pathway_to_receiver_col[pathway] = col
                    else:
                        pathway_to_receiver_col[col] = col
        
        # Combine all values and find global maximums for significant pathways only
        all_sender_df = pd.concat(all_sender_values)
        all_receiver_df = pd.concat(all_receiver_values)
        
        # Filter for significant pathways only
        significant_sender_cols = [pathway_to_sender_col[p] for p in significant_pathways if p in pathway_to_sender_col]
        significant_receiver_cols = [pathway_to_receiver_col[p] for p in significant_pathways if p in pathway_to_receiver_col]
        
        global_sender_max = all_sender_df[significant_sender_cols].max().max() if significant_sender_cols else 0
        global_receiver_max = all_receiver_df[significant_receiver_cols].max().max() if significant_receiver_cols else 0
        
        # Get valid pathways that have both sender and receiver mappings
        valid_pathways = [p for p in significant_pathways if p in pathway_to_sender_col or p in pathway_to_receiver_col]
        
        if not valid_pathways:
            logging.warning("No valid pathways found for plotting")
            return
        
        # Process each group
        for group_name, group_data in sample_categories.items():
            logging.info(f"Processing split plots for group: {group_name}")
            
            # Average the data across samples in the group
            avg_sender = pd.concat([d['sender'] for d in group_data.values()]).groupby(level=0).mean()
            avg_receiver = pd.concat([d['receiver'] for d in group_data.values()]).groupby(level=0).mean()
            
            # Get all celltypes
            celltypes = sorted(list(set(avg_sender.index) | set(avg_receiver.index)))
            
            # Filter pathways to those valid for this group
            group_valid_pathways = []
            for pathway in valid_pathways:
                sender_col = pathway_to_sender_col.get(pathway)
                receiver_col = pathway_to_receiver_col.get(pathway)
                
                if (sender_col and sender_col in avg_sender.columns) or \
                   (receiver_col and receiver_col in avg_receiver.columns):
                    group_valid_pathways.append(pathway)
            
            if not group_valid_pathways or not celltypes:
                continue
            
            # Create figure with two subplots side by side
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 10))
            
            # Create matrices for all pathways
            sender_matrix = np.zeros((len(group_valid_pathways), len(celltypes)))
            receiver_matrix = np.zeros((len(group_valid_pathways), len(celltypes)))
            
            # Create sender plot
            for i, pathway in enumerate(group_valid_pathways):
                sender_col = pathway_to_sender_col.get(pathway)
                if sender_col:
                    for j, celltype in enumerate(celltypes):
                        if celltype in avg_sender.index:
                            try:
                                val = avg_sender.loc[celltype, sender_col]
                                sender_matrix[i, j] = val / global_sender_max if global_sender_max > 0 else 0
                            except KeyError:
                                continue
            
            # Create receiver plot
            for i, pathway in enumerate(group_valid_pathways):
                receiver_col = pathway_to_receiver_col.get(pathway)
                if receiver_col:
                    for j, celltype in enumerate(celltypes):
                        if celltype in avg_receiver.index:
                            try:
                                val = avg_receiver.loc[celltype, receiver_col]
                                receiver_matrix[i, j] = val / global_receiver_max if global_receiver_max > 0 else 0
                            except KeyError:
                                continue
            
            # Plot sender heatmap
            im1 = ax1.imshow(sender_matrix, aspect='auto', cmap='viridis')
            ax1.set_xticks(range(len(celltypes)))
            ax1.set_yticks(range(len(group_valid_pathways)))
            ax1.set_xticklabels(celltypes, rotation=45, ha='right')
            ax1.set_yticklabels(group_valid_pathways)
            ax1.set_title(f'{group_name}\nSender Signals')
            plt.colorbar(im1, ax=ax1, label='Normalized Signal Strength')
            
            # Plot receiver heatmap
            im2 = ax2.imshow(receiver_matrix, aspect='auto', cmap='viridis')
            ax2.set_xticks(range(len(celltypes)))
            ax2.set_yticks(range(len(group_valid_pathways)))
            ax2.set_xticklabels(celltypes, rotation=45, ha='right')
            ax2.set_yticklabels(group_valid_pathways)
            ax2.set_title(f'{group_name}\nReceiver Signals')
            plt.colorbar(im2, ax=ax2, label='Normalized Signal Strength')
            
            # Adjust layout and save
            plt.tight_layout()
            safe_group = group_name.replace(" ", "_")
            filename = f'figures/signaling/{safe_group}_split_signals.png'
            fig.savefig(filename, bbox_inches='tight', dpi=300)
            fig.savefig(filename.replace('.png', '.svg'), bbox_inches='tight')
            logging.info(f"Saved split signal plot: {filename}")
            plt.close(fig)
            
    except Exception as e:
        logging.error(f"Error in create_split_signal_plots: {e}")
        logging.exception("Detailed traceback:")
        plt.close('all')

def create_patient_specific_comparison(results_by_celltype, patient_id, all_celltypes):
    """Create patient-specific comparison between biopsy and surgical samples.
    
    Args:
        results_by_celltype (dict): Dictionary containing results for all samples
        patient_id (str): ID of the patient to analyze
        all_celltypes (set): Set of all cell types to analyze
    """
    logging.info(f"Starting patient-specific comparison for {patient_id}")
    
    # Get biopsy and surgical samples for this patient
    biopsy_samples = {k: v for k, v in results_by_celltype.items() 
                     if patient_id in k and 'S1' in k}
    surgical_samples = {k: v for k, v in results_by_celltype.items() 
                       if patient_id in k and 'S2' in k}
    
    logging.info(f"Found {len(biopsy_samples)} biopsy and {len(surgical_samples)} surgical samples")
    print(f"\nDEBUG: Patient {patient_id}")
    print(f"Biopsy samples: {list(biopsy_samples.keys())}")
    print(f"Surgical samples: {list(surgical_samples.keys())}")
    print(f"Cell types to analyze: {sorted(all_celltypes)}")
    
    if not biopsy_samples or not surgical_samples:
        logging.warning(f"Could not find both biopsy and surgical samples for {patient_id}")
        return
    
    comparison_name = f"{patient_id} Biopsy vs Surgical"
    output_prefix = f"{patient_id.lower()}_biopsy_vs_surgical"
    
    patient_pathways_data = {}
    all_pathways = set()
    
    # Process each celltype
    for celltype in all_celltypes:
        logging.info(f"Analyzing celltype: {celltype}")
        
        pathway_results = create_comparison_plot(
            celltype,
            biopsy_samples,
            surgical_samples,
            comparison_name,
            output_prefix
        )
        
        if pathway_results:
            # Store all pathways, not just significant ones
            patient_pathways_data[celltype] = pathway_results
            all_pathways.update(pathway_results.keys())
            
            # Log pathway counts
            sig_pathways = sum(1 for p_data in pathway_results.values()
                             if min(p_data.get('sender', {}).get('p_val', 1.0),
                                  p_data.get('receiver', {}).get('p_val', 1.0)) < 0.05)
    
    logging.info(f"Total pathways found for {patient_id}: {len(all_pathways)}")
    
    # Create visualizations if we found any pathways
    if patient_pathways_data:
        create_pathway_bar_plots({comparison_name: patient_pathways_data}, plot_all=True)
        
        # Create split signal plots for all pathways
        if all_pathways:
            create_split_signal_plots(
                {
                    f"{patient_id} Biopsy": biopsy_samples,
                    f"{patient_id} Surgical": surgical_samples
                },
                all_pathways
            )
        else:
            logging.warning(f"No pathways found for {patient_id}")

def create_pathway_bar_plots(all_pathways_data, plot_all=False):
    """Create bar plots for each pathway comparison.
    
    Args:
        all_pathways_data (dict): Dictionary of pathway data by comparison
        plot_all (bool): If True, plot all pathways regardless of significance
    """
    try:
        # Close any existing figures at the start
        plt.close('all')
        
        # Process each comparison
        for comparison_name, celltype_data in all_pathways_data.items():
            logging.info(f"Creating bar plots for comparison: {comparison_name}")
            
            # Collect pathways data for this comparison
            pathways_data = []
            
            # Determine if this is a patient-specific comparison
            is_patient_comparison = "Biopsy vs Surgical" in comparison_name
            
            # New structure: celltype_data is a dictionary of celltypes to pathway dictionaries
            for celltype, pathways in celltype_data.items():
                for pathway_name, pathway_info in pathways.items():
                    # Get values for both sender and receiver
                    sender_data = pathway_info.get('sender', {})
                    receiver_data = pathway_info.get('receiver', {})
                    
                    if sender_data or receiver_data:
                        if sender_data:
                            pathways_data.append({
                                'pathway': pathway_name,
                                'celltype': celltype,
                                'type': 'sender',
                                'p_val': sender_data.get('p_val', 1.0),
                                'group1_val': sender_data.get('group1_mean', 0),
                                'group2_val': sender_data.get('group2_mean', 0),
                                'is_significant': sender_data.get('p_val', 1.0) < 0.05
                            })
                        
                        if receiver_data:
                            pathways_data.append({
                                'pathway': pathway_name,
                                'celltype': celltype,
                                'type': 'receiver',
                                'p_val': receiver_data.get('p_val', 1.0),
                                'group1_val': receiver_data.get('group1_mean', 0),
                                'group2_val': receiver_data.get('group2_mean', 0),
                                'is_significant': receiver_data.get('p_val', 1.0) < 0.05
                            })
            
            if not pathways_data:
                logging.info(f"No pathways found for {comparison_name}")
                continue
                
            # Sort by p-value and pathway name to group same pathways together
            pathways_data.sort(key=lambda x: (x['p_val'], x['pathway']))
            
            # Create figure with adjusted margins
            n_pathways = len(pathways_data)
            fig_height = max(8, n_pathways * 0.4)  # Adjust height based on number of pathways
            fig = plt.figure(figsize=(12, fig_height))
            
            # Create axes with adjusted position to accommodate labels
            ax = fig.add_axes([0.3, 0.1, 0.6, 0.85])  # [left, bottom, width, height]
            
            # Prepare data for plotting
            pathways = [f"{d['pathway']}\n({d['celltype']}, {d['type']})" for d in pathways_data]
            group1_vals = [d['group1_val'] for d in pathways_data]
            group2_vals = [d['group2_val'] for d in pathways_data]
            is_significant = [d['is_significant'] for d in pathways_data]
            types = [d['type'] for d in pathways_data]
            
            # Create bar positions
            y_pos = np.arange(len(pathways))
            bar_width = 0.35
            
            # Create bars with colors based on sender/receiver for patient comparisons
            # or significance status for other comparisons
            for i, (g1_val, g2_val, sig, type_) in enumerate(zip(group1_vals, group2_vals, is_significant, types)):
                if is_patient_comparison:
                    # For patient comparisons, use sender/receiver coloring
                    if type_ == 'sender':
                        color1 = 'skyblue'
                        color2 = 'lightblue'
                    else:  # receiver
                        color1 = 'lightcoral'
                        color2 = 'mistyrose'
                    alpha = 0.8
                else:
                    # For other comparisons, use significance coloring
                    color1 = 'skyblue' if sig else 'lightblue'
                    color2 = 'lightcoral' if sig else 'mistyrose'
                    alpha = 0.8 if sig else 0.6
                
                ax.barh(y_pos[i] - bar_width/2, g1_val, bar_width, color=color1, alpha=alpha)
                ax.barh(y_pos[i] + bar_width/2, g2_val, bar_width, color=color2, alpha=alpha)
            
            # Add appropriate legend based on comparison type
            if is_patient_comparison:
                legend_elements = [
                    plt.Rectangle((0,0),1,1, color='skyblue', alpha=0.8, label='Biopsy (Sender)'),
                    plt.Rectangle((0,0),1,1, color='lightblue', alpha=0.8, label='Surgical (Sender)'),
                    plt.Rectangle((0,0),1,1, color='lightcoral', alpha=0.8, label='Biopsy (Receiver)'),
                    plt.Rectangle((0,0),1,1, color='mistyrose', alpha=0.8, label='Surgical (Receiver)')
                ]
            else:
                legend_elements = [
                    plt.Rectangle((0,0),1,1, color='skyblue', alpha=0.8, label='Group 1 (Significant)'),
                    plt.Rectangle((0,0),1,1, color='lightblue', alpha=0.6, label='Group 1 (Not Significant)'),
                    plt.Rectangle((0,0),1,1, color='lightcoral', alpha=0.8, label='Group 2 (Significant)'),
                    plt.Rectangle((0,0),1,1, color='mistyrose', alpha=0.6, label='Group 2 (Not Significant)')
                ]
            ax.legend(handles=legend_elements, bbox_to_anchor=(1.02, 1), loc='upper left')
            
            # Customize plot
            ax.set_yticks(y_pos)
            ax.set_yticklabels(pathways)
            ax.invert_yaxis()  # Invert y-axis to show most significant at top
            ax.set_xlabel('Signal Strength')
            ax.set_title(comparison_name)
            
            # Add p-values
            max_val = max(max(group1_vals), max(group2_vals))
            for i, data in enumerate(pathways_data):
                text_color = 'black' if data['is_significant'] else 'gray'
                ax.text(max_val * 1.05, 
                       y_pos[i], 
                       f"p={data['p_val']:.2e}", 
                       va='center',
                       color=text_color)
            
            # Extend x-axis to accommodate p-values
            ax.set_xlim(0, max_val * 1.3)
            
            # Save plot
            safe_comparison = comparison_name.replace(" ", "_").replace(".", "").replace(",", "")
            filename = f'figures/signaling/{safe_comparison}_barplot.png'
            fig.savefig(filename, bbox_inches='tight', dpi=300)
            fig.savefig(filename.replace('.png', '.svg'), bbox_inches='tight')
            logging.info(f"Saved bar plot: {filename}")
            plt.close(fig)
            
    except Exception as e:
        logging.error(f"Error in create_pathway_bar_plots: {e}")
        logging.exception("Detailed traceback:")
        plt.close('all')

def plot_single_group_pathway(commot_results_all_samples, group_identifier_keyword, group_display_name, top_n=10, output_dir='figures/signaling/single_group_plots'):
    """
    Plot top N pathways for a single group prioritizing "unipolar" signals.
    
    Unipolar signals are pathways that have high activity in either sender OR receiver
    but not both - this captures more specific signaling patterns.
    
    Args:
        commot_results_all_samples (dict): Dictionary of COMMOT results by sample
        group_identifier_keyword (str): Keyword to identify samples (e.g., "S1" for biopsy)
        group_display_name (str): Display name for the group (e.g., "Biopsy Samples")
        top_n (int): Number of top pathways to plot per cell type
        output_dir (str): Directory to save plots
    
    Returns:
        dict: Dictionary with top pathways by cell type for sender and receiver
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import os
    import logging
    
    def extract_celltype(index_name):
        """Extract celltype name from index."""
        if '_genes_pass1' in index_name:
            index_name = index_name.replace('_genes_pass1', '')
        if '_' in index_name:
            index_name = index_name.split('_', 1)[1]
        return index_name
    
    # Filter samples by group identifier
    group_samples = {k: v for k, v in commot_results_all_samples.items() 
                     if group_identifier_keyword in k}
    
    if not group_samples:
        logging.warning(f"No samples found with identifier '{group_identifier_keyword}'")
        return None
    
    logging.info(f"Found {len(group_samples)} samples for {group_display_name}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Collect all pathways and cell types
    all_pathways = set()
    all_celltypes = set()
    
    for sample_data in group_samples.values():
        # Get pathways (remove prefixes)
        sender_pathways = {col.replace('s-', '') for col in sample_data['sender_signal'].columns}
        receiver_pathways = {col.replace('r-', '') for col in sample_data['receiver_signal'].columns}
        all_pathways.update(sender_pathways | receiver_pathways)
        
        # Get cell types
        sender_celltypes = set(sample_data['sender_signal'].index.map(extract_celltype))
        receiver_celltypes = set(sample_data['receiver_signal'].index.map(extract_celltype))
        all_celltypes.update(sender_celltypes | receiver_celltypes)
    
    # Remove total-total if present
    if 'total-total' in all_pathways:
        all_pathways.remove('total-total')
    
    logging.info(f"Found {len(all_pathways)} pathways and {len(all_celltypes)} cell types")
    
    # Storage for results
    results = {
        'top_sender_pathways_by_celltype': {},
        'top_receiver_pathways_by_celltype': {}
    }
    
    # Process each cell type
    for celltype in sorted(all_celltypes):
        logging.info(f"Processing cell type: {celltype}")
        
        pathway_scores = {}
        
        # Calculate unipolar scores for each pathway
        for pathway in all_pathways:
            sender_values = []
            receiver_values = []
            
            # Collect values across all samples
            for sample_data in group_samples.values():
                # Get sender values for this celltype and pathway
                sender_mask = sample_data['sender_signal'].index.map(extract_celltype) == celltype
                if any(sender_mask):
                    sender_col = f"s-{pathway}"
                    if sender_col in sample_data['sender_signal'].columns:
                        vals = sample_data['sender_signal'].loc[sender_mask, sender_col].values
                        sender_values.extend(vals)
                
                # Get receiver values for this celltype and pathway
                receiver_mask = sample_data['receiver_signal'].index.map(extract_celltype) == celltype
                if any(receiver_mask):
                    receiver_col = f"r-{pathway}"
                    if receiver_col in sample_data['receiver_signal'].columns:
                        vals = sample_data['receiver_signal'].loc[receiver_mask, receiver_col].values
                        receiver_values.extend(vals)
            
            if sender_values or receiver_values:
                # Calculate mean signals
                sender_mean = np.mean(sender_values) if sender_values else 0.0
                receiver_mean = np.mean(receiver_values) if receiver_values else 0.0
                
                # Calculate unipolar score: prioritize pathways with high activity in one direction
                # but low activity in the other direction
                # Score = max(sender, receiver) - min(sender, receiver)
                # This gives higher scores to pathways that are "unipolar"
                max_signal = max(sender_mean, receiver_mean)
                min_signal = min(sender_mean, receiver_mean)
                unipolar_score = max_signal - min_signal
                
                # Also consider the absolute magnitude to avoid tiny differences
                magnitude_weight = max_signal
                
                # Final score combines unipolar tendency with overall magnitude
                final_score = unipolar_score * magnitude_weight
                
                pathway_scores[pathway] = {
                    'sender_mean': sender_mean,
                    'receiver_mean': receiver_mean,
                    'unipolar_score': unipolar_score,
                    'magnitude_weight': magnitude_weight,
                    'final_score': final_score,
                    'dominant_type': 'sender' if sender_mean > receiver_mean else 'receiver'
                }
        
        # Sort pathways by unipolar score (descending)
        sorted_pathways = sorted(pathway_scores.items(), 
                               key=lambda x: x[1]['final_score'], 
                               reverse=True)
        
        # Separate into sender-dominant and receiver-dominant
        sender_dominant = [(p, s) for p, s in sorted_pathways if s['dominant_type'] == 'sender']
        receiver_dominant = [(p, s) for p, s in sorted_pathways if s['dominant_type'] == 'receiver']
        
        # Get top N for each type
        top_sender_pathways = [p for p, s in sender_dominant[:top_n]]
        top_receiver_pathways = [p for p, s in receiver_dominant[:top_n]]
        
        results['top_sender_pathways_by_celltype'][celltype] = top_sender_pathways
        results['top_receiver_pathways_by_celltype'][celltype] = top_receiver_pathways
        
        # Create plots for this cell type if there are pathways to plot
        if top_sender_pathways or top_receiver_pathways:
            # Create figure with subplots for sender and receiver
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
            
            # Plot top sender pathways
            if top_sender_pathways:
                sender_data = [(p, pathway_scores[p]) for p in top_sender_pathways]
                pathways_s = [p for p, _ in sender_data]
                scores_s = [s['sender_mean'] for _, s in sender_data]
                unipolar_s = [s['unipolar_score'] for _, s in sender_data]
                
                bars1 = ax1.barh(range(len(pathways_s)), scores_s, color='skyblue', alpha=0.8)
                ax1.set_yticks(range(len(pathways_s)))
                ax1.set_yticklabels(pathways_s)
                ax1.set_xlabel('Mean Sender Signal')
                ax1.set_title(f'{celltype} - Top {top_n} Unipolar Sender Pathways\n{group_display_name}')
                ax1.grid(True, alpha=0.3)
                
                # Set x-axis limits to accommodate the U text
                max_score = max(scores_s) if scores_s else 1
                ax1.set_xlim(0, max_score * 1.2)
                
                # Add unipolar scores as text
                for i, (bar, unipolar_score) in enumerate(zip(bars1, unipolar_s)):
                    ax1.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                           f'U:{unipolar_score:.3f}', va='center', fontsize=14)
            else:
                ax1.text(0.5, 0.5, 'No sender-dominant pathways', 
                        ha='center', va='center', transform=ax1.transAxes)
                ax1.set_title(f'{celltype} - Top {top_n} Unipolar Sender Pathways\n{group_display_name}')
            
            # Plot top receiver pathways
            if top_receiver_pathways:
                receiver_data = [(p, pathway_scores[p]) for p in top_receiver_pathways]
                pathways_r = [p for p, _ in receiver_data]
                scores_r = [s['receiver_mean'] for _, s in receiver_data]
                unipolar_r = [s['unipolar_score'] for _, s in receiver_data]
                
                bars2 = ax2.barh(range(len(pathways_r)), scores_r, color='lightcoral', alpha=0.8)
                ax2.set_yticks(range(len(pathways_r)))
                ax2.set_yticklabels(pathways_r)
                ax2.set_xlabel('Mean Receiver Signal')
                ax2.set_title(f'{celltype} - Top {top_n} Unipolar Receiver Pathways\n{group_display_name}')
                ax2.grid(True, alpha=0.3)
                
                # Set x-axis limits to accommodate the U text
                max_score = max(scores_r) if scores_r else 1
                ax2.set_xlim(0, max_score * 1.2)
                
                # Add unipolar scores as text
                for i, (bar, unipolar_score) in enumerate(zip(bars2, unipolar_r)):
                    ax2.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                           f'U:{unipolar_score:.3f}', va='center', fontsize=14)
            else:
                ax2.text(0.5, 0.5, 'No receiver-dominant pathways', 
                        ha='center', va='center', transform=ax2.transAxes)
                ax2.set_title(f'{celltype} - Top {top_n} Unipolar Receiver Pathways\n{group_display_name}')
            
            # Add explanation text
            fig.text(0.02, 0.02, 'U: Unipolar Score (higher = more unipolar)', fontsize=14)
            
            plt.tight_layout()
            
            # Save plot
            safe_celltype = celltype.replace('/', '_').replace(' ', '_')
            filename = os.path.join(output_dir, f'{safe_celltype}_unipolar_pathways.png')
            fig.savefig(filename, bbox_inches='tight', dpi=300)
            fig.savefig(filename.replace('.png', '.svg'), bbox_inches='tight')
            plt.close(fig)
            
            logging.info(f"Saved unipolar pathway plot for {celltype}: {filename}")
        
        # Log top pathways for this cell type
        if top_sender_pathways:
            logging.info(f"Top sender pathways for {celltype}: {top_sender_pathways[:3]}...")
        if top_receiver_pathways:
            logging.info(f"Top receiver pathways for {celltype}: {top_receiver_pathways[:3]}...")
    
    # Create summary plot showing unipolar scores across all cell types
    if results['top_sender_pathways_by_celltype'] or results['top_receiver_pathways_by_celltype']:
        create_unipolar_summary_plot(pathway_scores, all_celltypes, group_display_name, output_dir)
    
    return results


def create_unipolar_summary_plot(all_pathway_scores, celltypes, group_display_name, output_dir):
    """Create a summary plot showing the most unipolar pathways across all cell types."""
    import matplotlib.pyplot as plt
    import numpy as np
    
    # This would need to be implemented if you want a cross-celltype summary
    # For now, just log that individual plots were created
    logging.info(f"Created individual unipolar pathway plots for {len(celltypes)} cell types in {group_display_name}")


def plot_filtered_pathway_analysis(commot_results_all_samples, celltype_filter_keywords, pathway_filter_keywords, 
                                  group_display_name, top_n=None, output_dir='figures/signaling/filtered_analysis'):
    """
    Plot pathways for specific cell types and pathway patterns (e.g., macrophages with IL/TGF pathways).
    Plots ALL found pathways sorted by highest magnitude.
    
    Args:
        commot_results_all_samples (dict): Dictionary of COMMOT results by sample
        celltype_filter_keywords (list): Keywords to filter cell types (e.g., ['macrophage'])
        pathway_filter_keywords (list): Keywords to filter pathways (e.g., ['IL', 'TGF'])
        group_display_name (str): Display name for the group (e.g., "Macrophage IL/TGF Analysis")
        top_n (int): Ignored - plots all pathways found (kept for compatibility)
        output_dir (str): Directory to save plots
    
    Returns:
        dict: Dictionary with all pathways by cell type for sender and receiver
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import os
    import logging
    
    def extract_celltype(index_name):
        """Extract celltype name from index."""
        if '_genes_pass1' in index_name:
            index_name = index_name.replace('_genes_pass1', '')
        if '_' in index_name:
            index_name = index_name.split('_', 1)[1]
        return index_name
    
    def matches_filter(text, filter_keywords):
        """Check if text contains any of the filter keywords (case-insensitive)."""
        text_lower = text.lower()
        return any(keyword.lower() in text_lower for keyword in filter_keywords)
    
    logging.info(f"Starting filtered pathway analysis for {group_display_name}")
    logging.info(f"Cell type filters: {celltype_filter_keywords}")
    logging.info(f"Pathway filters: {pathway_filter_keywords}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Collect all pathways and cell types that match filters
    filtered_pathways = set()
    filtered_celltypes = set()
    
    for sample_data in commot_results_all_samples.values():
        # Get pathways that match filter (remove prefixes)
        sender_pathways = {col.replace('s-', '') for col in sample_data['sender_signal'].columns
                          if matches_filter(col.replace('s-', ''), pathway_filter_keywords)}
        receiver_pathways = {col.replace('r-', '') for col in sample_data['receiver_signal'].columns
                           if matches_filter(col.replace('r-', ''), pathway_filter_keywords)}
        filtered_pathways.update(sender_pathways | receiver_pathways)
        
        # Get cell types that match filter
        sender_celltypes = {extract_celltype(idx) for idx in sample_data['sender_signal'].index
                           if matches_filter(extract_celltype(idx), celltype_filter_keywords)}
        receiver_celltypes = {extract_celltype(idx) for idx in sample_data['receiver_signal'].index
                            if matches_filter(extract_celltype(idx), celltype_filter_keywords)}
        filtered_celltypes.update(sender_celltypes | receiver_celltypes)
    
    # Remove total-total if present
    if 'total-total' in filtered_pathways:
        filtered_pathways.remove('total-total')
    
    logging.info(f"Found {len(filtered_pathways)} matching pathways: {sorted(filtered_pathways)}")
    logging.info(f"Found {len(filtered_celltypes)} matching cell types: {sorted(filtered_celltypes)}")
    
    if not filtered_pathways or not filtered_celltypes:
        logging.warning("No matching pathways or cell types found with the given filters")
        return None
    
    # Storage for results
    results = {
        'all_sender_pathways_by_celltype': {},
        'all_receiver_pathways_by_celltype': {},
        'filtered_pathways': list(filtered_pathways),
        'filtered_celltypes': list(filtered_celltypes)
    }
    
    # Process each filtered cell type
    for celltype in sorted(filtered_celltypes):
        logging.info(f"Processing cell type: {celltype}")
        
        pathway_scores = {}
        
        # Calculate scores for each filtered pathway
        for pathway in filtered_pathways:
            sender_values = []
            receiver_values = []
            
            # Collect values across all samples
            for sample_data in commot_results_all_samples.values():
                # Get sender values for this celltype and pathway
                sender_mask = sample_data['sender_signal'].index.map(extract_celltype) == celltype
                if any(sender_mask):
                    sender_col = f"s-{pathway}"
                    if sender_col in sample_data['sender_signal'].columns:
                        vals = sample_data['sender_signal'].loc[sender_mask, sender_col].values
                        sender_values.extend(vals)
                
                # Get receiver values for this celltype and pathway
                receiver_mask = sample_data['receiver_signal'].index.map(extract_celltype) == celltype
                if any(receiver_mask):
                    receiver_col = f"r-{pathway}"
                    if receiver_col in sample_data['receiver_signal'].columns:
                        vals = sample_data['receiver_signal'].loc[receiver_mask, receiver_col].values
                        receiver_values.extend(vals)
            
            if sender_values or receiver_values:
                # Calculate mean signals
                sender_mean = np.mean(sender_values) if sender_values else 0.0
                receiver_mean = np.mean(receiver_values) if receiver_values else 0.0
                
                # Calculate magnitude (highest signal strength)
                max_signal = max(sender_mean, receiver_mean)
                
                pathway_scores[pathway] = {
                    'sender_mean': sender_mean,
                    'receiver_mean': receiver_mean,
                    'magnitude': max_signal,
                    'dominant_type': 'sender' if sender_mean > receiver_mean else 'receiver'
                }
        
        if not pathway_scores:
            logging.warning(f"No pathway data found for cell type {celltype}")
            continue
        
        # Sort pathways by magnitude (descending) - highest signal strength first
        sorted_pathways = sorted(pathway_scores.items(), 
                               key=lambda x: x[1]['magnitude'], 
                               reverse=True)
        
        # Separate into sender-dominant and receiver-dominant (but keep ALL, not just top N)
        sender_dominant = [(p, s) for p, s in sorted_pathways if s['dominant_type'] == 'sender']
        receiver_dominant = [(p, s) for p, s in sorted_pathways if s['dominant_type'] == 'receiver']
        
        # Get ALL pathways for each type (no filtering by top_n)
        all_sender_pathways = [p for p, s in sender_dominant]
        all_receiver_pathways = [p for p, s in receiver_dominant]
        
        results['all_sender_pathways_by_celltype'][celltype] = all_sender_pathways
        results['all_receiver_pathways_by_celltype'][celltype] = all_receiver_pathways
        
        # Create plots for this cell type if there are pathways to plot
        if all_sender_pathways or all_receiver_pathways:
            # Create figure with subplots for sender and receiver
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
            
            # Plot all sender pathways
            if all_sender_pathways:
                sender_data = [(p, pathway_scores[p]) for p in all_sender_pathways]
                pathways_s = [p for p, _ in sender_data]
                scores_s = [s['sender_mean'] for _, s in sender_data]
                magnitudes_s = [s['magnitude'] for _, s in sender_data]
                
                bars1 = ax1.barh(range(len(pathways_s)), scores_s, color='skyblue', alpha=0.8)
                ax1.set_yticks(range(len(pathways_s)))
                ax1.set_yticklabels(pathways_s)
                ax1.set_xlabel('Mean Sender Signal')
                ax1.set_title(f'{celltype} - All IL/TGF Sender Pathways\n{group_display_name}')
                ax1.grid(True, alpha=0.3)
                
                # Set x-axis limits to accommodate the magnitude text
                max_score = max(scores_s) if scores_s else 1
                ax1.set_xlim(0, max_score * 1.2)
                
                # Add magnitude scores as text
                for i, (bar, magnitude) in enumerate(zip(bars1, magnitudes_s)):
                    ax1.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                           f'M:{magnitude:.3f}', va='center', fontsize=14)
            else:
                ax1.text(0.5, 0.5, 'No sender-dominant pathways', 
                        ha='center', va='center', transform=ax1.transAxes)
                ax1.set_title(f'{celltype} - All IL/TGF Sender Pathways\n{group_display_name}')
            
            # Plot all receiver pathways
            if all_receiver_pathways:
                receiver_data = [(p, pathway_scores[p]) for p in all_receiver_pathways]
                pathways_r = [p for p, _ in receiver_data]
                scores_r = [s['receiver_mean'] for _, s in receiver_data]
                magnitudes_r = [s['magnitude'] for _, s in receiver_data]
                
                bars2 = ax2.barh(range(len(pathways_r)), scores_r, color='lightcoral', alpha=0.8)
                ax2.set_yticks(range(len(pathways_r)))
                ax2.set_yticklabels(pathways_r)
                ax2.set_xlabel('Mean Receiver Signal')
                ax2.set_title(f'{celltype} - All IL/TGF Receiver Pathways\n{group_display_name}')
                ax2.grid(True, alpha=0.3)
                
                # Set x-axis limits to accommodate the magnitude text
                max_score = max(scores_r) if scores_r else 1
                ax2.set_xlim(0, max_score * 1.2)
                
                # Add magnitude scores as text
                for i, (bar, magnitude) in enumerate(zip(bars2, magnitudes_r)):
                    ax2.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                           f'M:{magnitude:.3f}', va='center', fontsize=14)
            else:
                ax2.text(0.5, 0.5, 'No receiver-dominant pathways', 
                        ha='center', va='center', transform=ax2.transAxes)
                ax2.set_title(f'{celltype} - All IL/TGF Receiver Pathways\n{group_display_name}')
            
            # Add explanation text
            fig.text(0.02, 0.02, 'M: Magnitude Score (signal strength)', fontsize=14)
            
            plt.tight_layout()
            
            # Save plot
            safe_celltype = celltype.replace('/', '_').replace(' ', '_')
            filename = os.path.join(output_dir, f'{safe_celltype}_IL_TGF_pathways.png')
            fig.savefig(filename, bbox_inches='tight', dpi=300)
            fig.savefig(filename.replace('.png', '.svg'), bbox_inches='tight')
            plt.close(fig)
            
            logging.info(f"Saved IL/TGF pathway plot for {celltype}: {filename}")
        
        # Log all pathways for this cell type
        if all_sender_pathways:
            logging.info(f"All IL/TGF sender pathways for {celltype}: {all_sender_pathways}")
        if all_receiver_pathways:
            logging.info(f"All IL/TGF receiver pathways for {celltype}: {all_receiver_pathways}")
    
    return results