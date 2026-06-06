#!/usr/bin/env python3
# PROVENANCE: Ported 2026-06-05 from
#   /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/05_commot_analysis.py
# Drives COMMOT cell-cell signaling analysis on the spatial biopsy samples
# and renders the per-cell-type sender_top_pathways SVGs used as Fig 7I.
# Calls commot_lib/utils.plot_single_group_pathway -> produces
#   figures/signaling/neil_aging_paper_summary/<CellType>_sender_top_pathways.{svg,png}
import os
import sys
import gc
import pickle
import json
import logging
import psutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from citegeist.core.analysis_functions import expand_prop_gex_adata
from utils import process_commot_analysis, standardize_adata, standardize_cell_profiles, standardize_layer_name, log_memory_usage, find_significant_pathways, create_comparison_plot, create_pathway_bar_plots, create_split_signal_plots, create_patient_specific_comparison, plot_single_group_pathway, plot_filtered_pathway_analysis   

# Create necessary directories first
os.makedirs('logs', exist_ok=True)
os.makedirs('data', exist_ok=True)
os.makedirs('figures', exist_ok=True)
os.makedirs('figures/signaling', exist_ok=True)

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/05_commot_analysis.log'),
        logging.StreamHandler(sys.stderr)
    ]
)


# Force stdout to flush immediately
sys.stdout.reconfigure(line_buffering=True)

logging.info("Starting COMMOT analysis script")
log_memory_usage()

# Load saved data
logging.info("Loading processed data...")
with open('data/biopsy_adatas.pkl', 'rb') as f:
    biopsy_adatas = pickle.load(f)
with open('data/surg_adatas.pkl', 'rb') as f:
    surg_adatas = pickle.load(f)
with open('data/cell_profiles.json', 'r') as f:
    cell_profiles = json.load(f)

RUN = False


if RUN:
    # Standardize cell profile keys to match the standardized layer names
    standardized_cell_profiles = standardize_cell_profiles(cell_profiles)

    # Create a mapping between standardized and original cell type names
    reverse_name_mapping = {}
    for original_name in cell_profiles.keys():
        std_name = standardize_layer_name(original_name)
        reverse_name_mapping[std_name] = original_name

    # COMMOT analysis using our utility function - process one sample at a time
    logging.info("Running COMMOT analysis...")
    commot_results = {}
    all_adata_dicts = [biopsy_adatas, surg_adatas]

    for adata_dict in all_adata_dicts:
        for sample_name, adata in adata_dict.items():
            logging.info(f"Processing COMMOT for {sample_name}...")
            try:
                # Standardize the AnnData object
                adata_copy = adata.copy()
                adata_copy, _ = standardize_adata(adata_copy, cell_profiles)
                
                # Expand data with standardized names
                expanded_adata = expand_prop_gex_adata(adata_copy, celltype_profile_dict=standardized_cell_profiles)
                
                # Clear copy
                del adata_copy
                gc.collect()
                
                commot_results[sample_name] = process_commot_analysis(expanded_adata, standardized_cell_profiles)
                del expanded_adata
                gc.collect()
                log_memory_usage()
            except Exception as e:
                logging.error(f"Error running COMMOT on {sample_name}: {e}")
                continue

    # Process COMMOT results
    logging.info("Processing COMMOT results...")
    logging.info("Commot results: %s", commot_results)
    # Save processed COMMOT results
    logging.info("Saving COMMOT results...")
    with open('data/commot_results.pkl', 'wb') as f:
        pickle.dump(commot_results, f)

# Load the saved results for analysis
logging.info("Loading saved COMMOT results for cross-comparison analysis...")
with open('data/commot_results.pkl', 'rb') as f:
    commot_results = pickle.load(f)

from utils import find_significant_pathways_by_comparison, plot_significant_results_by_celltype

# Define comparisons list with proper format
comparisons = [
    ('non_progressing_biopsy', 'non_progressing_surgical', 'Non-progressing Biopsy vs Surgical', 'np_biopsy_vs_surgical'),
    ('progressing_biopsy', 'progressing_surgical', 'Progressing Biopsy vs Surgical', 'p_biopsy_vs_surgical'),
    ('non_progressing_biopsy', 'progressing_biopsy', 'Non-progressing vs Progressing Biopsy', 'np_vs_p_biopsy'),
    ('non_progressing_surgical', 'progressing_surgical', 'Non-progressing vs Progressing Surgical', 'np_vs_p_surgical')
]

# Separate samples into different categories for analysis
sample_categories = {
    'non_progressing_biopsy': {k: v for k, v in commot_results.items() 
                              if not any(x in k for x in ['P1', 'P4']) and 'S1' in k},
    'non_progressing_surgical': {k: v for k, v in commot_results.items() 
                                if not any(x in k for x in ['P1', 'P4']) and 'S2' in k},
    'progressing_biopsy': {k: v for k, v in commot_results.items() 
                          if any(x in k for x in ['P1', 'P4']) and 'S1' in k},
    'progressing_surgical': {k: v for k, v in commot_results.items() 
                           if any(x in k for x in ['P1', 'P4']) and 'S2' in k}
}

print(sample_categories)



for comparison in comparisons:
    significant_results = find_significant_pathways_by_comparison(sample_categories, comparison)
    plot_significant_results_by_celltype(
        significant_results, 
        comparison,
        output_dir='figures/signaling/celltype_comparisons'
    )
    print(significant_results)


# For Neil's aging paper, gather together all the biopsy sample commot results and plot the top 10 pathways by average signal strength, one for reciever and one for sender. 
neil_aging_plot_dir = 'figures/signaling/neil_aging_paper_summary'
os.makedirs(neil_aging_plot_dir, exist_ok=True)

# Call the new function for biopsy samples (samples with "S1" identifier)
logging.info("Generating plots for Neil's aging paper...")
biopsy_results = plot_single_group_pathway(
    commot_results_all_samples=commot_results,
    group_identifier_keyword="S1",  # "S1" indicates biopsy samples
    group_display_name="Biopsy Samples",
    top_n=10,
    output_dir=neil_aging_plot_dir
)

if biopsy_results:
    # Log a summary of the top pathways for each cell type
    sender_results = biopsy_results.get('top_sender_pathways_by_celltype', {})
    receiver_results = biopsy_results.get('top_receiver_pathways_by_celltype', {})
    
    logging.info(f"Generated per-cell type plots for {len(sender_results)} cell types (sender signals)")
    logging.info(f"Generated per-cell type plots for {len(receiver_results)} cell types (receiver signals)")
    
    # Log details for a few cell types as examples
    for celltype, pathways in list(sender_results.items())[:3]:  # Show first 3 cell types
        if pathways:
            logging.info(f"Top sender pathways for {celltype}: {pathways[:3]}...")
    
    for celltype, pathways in list(receiver_results.items())[:3]:  # Show first 3 cell types
        if pathways:
            logging.info(f"Top receiver pathways for {celltype}: {pathways[:3]}...")
else:
    logging.warning("No results were returned from the plot_single_group_pathway function")

# Macrophage IL/TGF pathway analysis - filter for macrophage cell types and IL/TGF pathways
macrophage_il_tgf_plot_dir = 'figures/signaling/macrophage_il_tgf_analysis'
os.makedirs(macrophage_il_tgf_plot_dir, exist_ok=True)

logging.info("Generating macrophage IL/TGF pathway analysis...")
macrophage_il_tgf_results = plot_filtered_pathway_analysis(
    commot_results_all_samples=commot_results,
    celltype_filter_keywords=["macrophage"],  # Filter for any celltype containing "macrophage"
    pathway_filter_keywords=["IL", "TGF"],    # Filter for pathways containing "IL" or "TGF"
    group_display_name="Macrophage IL/TGF Signaling",
    output_dir=macrophage_il_tgf_plot_dir
)

if macrophage_il_tgf_results:
    # Log a summary of the results
    sender_results = macrophage_il_tgf_results.get('all_sender_pathways_by_celltype', {})
    receiver_results = macrophage_il_tgf_results.get('all_receiver_pathways_by_celltype', {})
    filtered_pathways = macrophage_il_tgf_results.get('filtered_pathways', [])
    filtered_celltypes = macrophage_il_tgf_results.get('filtered_celltypes', [])
    
    logging.info(f"Macrophage IL/TGF Analysis Results:")
    logging.info(f"- Found {len(filtered_celltypes)} macrophage cell types: {filtered_celltypes}")
    logging.info(f"- Found {len(filtered_pathways)} IL/TGF pathways: {filtered_pathways}")
    logging.info(f"- Generated plots for {len(sender_results)} cell types (sender signals)")
    logging.info(f"- Generated plots for {len(receiver_results)} cell types (receiver signals)")
    
    # Log all pathways for each macrophage cell type
    for celltype in filtered_celltypes:
        sender_pathways = sender_results.get(celltype, [])
        receiver_pathways = receiver_results.get(celltype, [])
        if sender_pathways:
            logging.info(f"All IL/TGF sender pathways for {celltype}: {sender_pathways}")
        if receiver_pathways:
            logging.info(f"All IL/TGF receiver pathways for {celltype}: {receiver_pathways}")
else:
    logging.warning("No results were returned from the macrophage IL/TGF pathway analysis")

sys.exit(0)

# Continue with the rest of the script
results_by_celltype = {}
for sample_name, results in commot_results.items():
    try:
        # Extract celltype names by removing both the barcode prefix and '_genes_pass1' suffix
        def extract_celltype(index_name):
            # Remove '_genes_pass1' suffix if present
            if '_genes_pass1' in index_name:
                index_name = index_name.replace('_genes_pass1', '')
            # Remove barcode prefix (everything before and including first '_')
            if '_' in index_name:
                index_name = index_name.split('_', 1)[1]
            return index_name

        sender_celltype = results['sender_signal'].index.map(extract_celltype)
        receiver_celltype = results['receiver_signal'].index.map(extract_celltype)
        
        sender_by_celltype = results['sender_signal'].groupby(sender_celltype).mean()
        receiver_by_celltype = results['receiver_signal'].groupby(receiver_celltype).mean()
        
        results_by_celltype[sample_name] = {
            'sender': sender_by_celltype,
            'receiver': receiver_by_celltype
        }
        gc.collect()
        log_memory_usage()

        # Log the number of pathways for each sample
        sender_pathways = set(col.replace('s-', '') for col in results['sender_signal'].columns)
        receiver_pathways = set(col.replace('r-', '') for col in results['receiver_signal'].columns)
        logging.info(f"Sample {sample_name}:")
        logging.info(f"  Sender pathways: {len(sender_pathways)}")
        logging.info(f"  Receiver pathways: {len(receiver_pathways)}")
        logging.info(f"  Total unique pathways: {len(sender_pathways | receiver_pathways)}")
    except Exception as e:
        logging.error(f"Error processing COMMOT results for {sample_name}: {e}")
        continue

# Clear commot results to save memory
del commot_results
gc.collect()

print(results_by_celltype)

# Save processed COMMOT results
logging.info("Saving COMMOT results by celltype...")
with open('data/commot_results_by_celltype.pkl', 'wb') as f:
    pickle.dump(results_by_celltype, f)

logging.info("COMMOT analysis complete!")

# Load the saved results for analysis
logging.info("Loading saved COMMOT results for cross-comparison analysis...")
with open('data/commot_results_by_celltype.pkl', 'rb') as f:
    results_by_celltype = pickle.load(f)

# Separate samples into different categories for analysis
sample_categories = {
    'non_progressing_biopsy': {k: v for k, v in results_by_celltype.items() 
                              if not any(x in k for x in ['P1', 'P4']) and 'S1' in k},
    'non_progressing_surgical': {k: v for k, v in results_by_celltype.items() 
                                if not any(x in k for x in ['P1', 'P4']) and 'S2' in k},
    'progressing_biopsy': {k: v for k, v in results_by_celltype.items() 
                          if any(x in k for x in ['P1', 'P4']) and 'S1' in k},
    'progressing_surgical': {k: v for k, v in results_by_celltype.items() 
                           if any(x in k for x in ['P1', 'P4']) and 'S2' in k}
}

# Get all unique celltypes from both the cell profiles and the results
all_celltypes = set(cell_profiles.keys())  # Get celltypes from original profiles only
logging.info(f"Using cell types from profiles: {sorted(all_celltypes)}")

logging.info(f"Sample categories: {sample_categories}")

# Create function that takes the comparison and then tests for significance on celltype by celltype level (in utils.py)
from utils import find_significant_pathways_by_comparison


for comparison in comparisons:
    find_significant_pathways_by_comparison(sample_categories, all_celltypes, comparison)




# First find all significant pathways by comparison
significant_pathways_by_comparison = find_significant_pathways(sample_categories, all_celltypes)
all_significant_pathways = set()
for pathways in significant_pathways_by_comparison.values():
    all_significant_pathways.update(pathways)

logging.info(f"Found significant pathways by comparison:")
for comparison, pathways in significant_pathways_by_comparison.items():
    logging.info(f"  {comparison}: {len(pathways)} pathways")
    logging.info(f"  Pathways: {sorted(pathways)}")
logging.info(f"Total unique significant pathways: {len(all_significant_pathways)}")

if not all_significant_pathways:
    logging.error("No significant pathways found. Check significance threshold and data.")
    sys.exit(1)

# Dictionary to store significant pathways found in bar plots for each comparison
bar_plot_significant_pathways = {}

# Initialize dictionary to store all pathway data
all_pathways_data = {}

# Define comparisons list with proper format
comparisons = [
    ('non_progressing_biopsy', 'non_progressing_surgical', 'Non-progressing Biopsy vs Surgical', 'np_biopsy_vs_surgical'),
    ('progressing_biopsy', 'progressing_surgical', 'Progressing Biopsy vs Surgical', 'p_biopsy_vs_surgical'),
    ('non_progressing_biopsy', 'progressing_biopsy', 'Non-progressing vs Progressing Biopsy', 'np_vs_p_biopsy'),
    ('non_progressing_surgical', 'progressing_surgical', 'Non-progressing vs Progressing Surgical', 'np_vs_p_surgical')
]

# Then create comparison bar plots
for data1_key, data2_key, comparison_name, output_prefix in comparisons:
    if data1_key in sample_categories and data2_key in sample_categories:
        logging.info(f"Processing comparison: {comparison_name}")
        comparison_key = f"{data1_key}_vs_{data2_key}"
        
        # Initialize storage for this comparison
        all_pathways_data[comparison_name] = {}
        bar_plot_significant_pathways[comparison_key] = set()
        
        # Process each celltype
        for celltype in all_celltypes:
            logging.info(f"Processing cell type: {celltype}")
            pathway_results = create_comparison_plot(
                celltype, 
                sample_categories[data1_key], 
                sample_categories[data2_key], 
                comparison_name, 
                output_prefix
            )
            
            # Store pathway results and collect significant pathways
            if pathway_results:
                all_pathways_data[comparison_name][celltype] = pathway_results
                # Add pathways to the set for this comparison
                for pathway_name, pathway_info in pathway_results.items():
                    sender_p_val = pathway_info.get('sender', {}).get('p_val', 1.0)
                    receiver_p_val = pathway_info.get('receiver', {}).get('p_val', 1.0)
                    if min(sender_p_val, receiver_p_val) < 0.05:
                        bar_plot_significant_pathways[comparison_key].add(pathway_name)

# Create bar plots for significant pathways
if all_pathways_data:
    create_pathway_bar_plots(all_pathways_data)
else:
    logging.error("No pathway data collected for bar plots.")

# Log the pathways found in bar plots
for comparison_key, pathways in bar_plot_significant_pathways.items():
    logging.info(f"Significant pathways found in bar plots for {comparison_key}: {len(pathways)}")
    logging.info(f"Pathways: {sorted(pathways)}")

# Create split signal plots using pathways found in bar plots
for data1_key, data2_key, comparison_name, output_prefix in comparisons:
    if data1_key in sample_categories and data2_key in sample_categories:
        comparison_key = f"{data1_key}_vs_{data2_key}"
        significant_pathways = bar_plot_significant_pathways[comparison_key]
        if significant_pathways:
            logging.info(f"Creating split signal plots for {comparison_key} with {len(significant_pathways)} pathways")
            create_split_signal_plots(
                {data1_key: sample_categories[data1_key], 
                 data2_key: sample_categories[data2_key]}, 
                significant_pathways
            )
        else:
            logging.warning(f"No significant pathways found in bar plots for comparison {comparison_name}")

logging.info("Running patient-specific comparisons...")
for patient_id in ['P1', 'P4']:
    create_patient_specific_comparison(results_by_celltype, patient_id, all_celltypes)

logging.info("Patient-specific analysis complete!")

logging.info("Analysis complete!")

