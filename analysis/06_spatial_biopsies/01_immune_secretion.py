import os
import sys
import pickle
import logging
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
from utils import log_memory_usage
from figure_config import setup_figure_params

# Create necessary directories
os.makedirs('logs', exist_ok=True)
os.makedirs('figures', exist_ok=True)
os.makedirs('figures/immune_secretion', exist_ok=True)

# Apply standardized figure parameters (14pt minimum, Arial font)
setup_figure_params()

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/01_immune_secretion.log'),
        logging.StreamHandler(sys.stderr)
    ]
)

# Force stdout to flush immediately
sys.stdout.reconfigure(line_buffering=True)

logging.info("Starting immune secretion analysis script")
log_memory_usage()

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

# Define cytokines of interest
cytokines = ['IL4', 'IL10', 'IL13', 'TGFB1', 'TGFB2', 'TGFB3']
logging.info(f"Looking for cytokines: {cytokines}")

# Macrophage-related keywords to identify relevant layers
macrophage_keywords = ['macrophage', 'Macrophage', 'CD163', 'CD14', 'HLA-DR', 'CD11c']

def find_macrophage_layers(adata):
    """Find layers that contain macrophage-related information"""
    macrophage_layers = []
    
    # Check layers
    if hasattr(adata, 'layers') and adata.layers:
        for layer_name in adata.layers.keys():
            if any(keyword in layer_name for keyword in macrophage_keywords):
                macrophage_layers.append(layer_name)
    
    # Check obs columns for cell type annotations
    macrophage_obs = []
    for col in adata.obs.columns:
        if any(keyword in col for keyword in macrophage_keywords):
            macrophage_obs.append(col)
    
    return macrophage_layers, macrophage_obs

def create_spatial_plot(adata, gene, layer_name=None, sample_name="", output_dir="figures/immune_secretion"):
    """Create and save spatial plot for a specific gene using scanpy"""
    try:
        # Check if gene exists
        if gene not in adata.var_names:
            logging.warning(f"Gene {gene} not found in var_names for {sample_name}")
            return False
        
        # Check if we have spatial coordinates
        if 'spatial' not in adata.obsm:
            logging.warning(f"No spatial coordinates found for {sample_name}")
            return False
        
        # Determine the plot title and filename
        if layer_name and layer_name in adata.layers:
            plot_title = f"{gene} expression in {layer_name} - {sample_name}"
            filename = f"{sample_name}_{gene}_{layer_name.replace(' ', '_')}_spatial"
            # Temporarily set the layer as the main expression for plotting
            adata_temp = adata.copy()
            adata_temp.X = adata_temp.layers[layer_name]
        else:
            plot_title = f"{gene} expression - {sample_name}"
            filename = f"{sample_name}_{gene}_spatial"
            adata_temp = adata
        
        # Create spatial plot using scanpy
        sc.pl.spatial(
            adata_temp,
            color=gene,
            title=plot_title,
            save=f"_{filename}.png",
            show=False,
            color_map='viridis'
        )
        
        # The file is saved by scanpy in the current working directory under 'figures/'
        # Move it to our desired output directory if different
        desired_path = os.path.join(output_dir, f"{filename}.png")

        # Try multiple possible scanpy output paths (prefix varies across scanpy versions)
        from pathlib import Path
        possible_paths = [
            Path(f"figures/show_{filename}.png"),
            Path(f"figures/{filename}.png"),
            Path(output_dir) / f"show_{filename}.png",
        ]
        saved_path = next((p for p in possible_paths if p.exists()), None)

        if saved_path is not None and str(saved_path) != desired_path:
            os.makedirs(output_dir, exist_ok=True)
            os.rename(str(saved_path), desired_path)
            logging.info(f"Saved spatial plot: {desired_path}")
        elif os.path.exists(desired_path):
            logging.info(f"Saved spatial plot: {desired_path}")
        else:
            logging.warning(f"Could not find saved plot file for {gene} in {sample_name}")
            return False

        # Save SVG alongside PNG
        svg_path = desired_path.replace('.png', '.svg')
        plt.savefig(svg_path, format='svg', bbox_inches='tight')
        logging.info(f"Saved SVG spatial plot: {svg_path}")

        return True
        
    except Exception as e:
        logging.error(f"Error creating spatial plot for {gene} in {sample_name}: {e}")
        return False

# Process each biopsy sample
for sample_name, adata in biopsy_adatas.items():
    logging.info(f"\nProcessing sample: {sample_name}")
    
    # Check which cytokines are present
    missing = [g for g in cytokines if g not in adata.var_names]
    present = [g for g in cytokines if g in adata.var_names]
    
    logging.info(f"Present cytokines: {present}")
    logging.info(f"Missing cytokines: {missing}")
    
    if not present:
        logging.warning(f"No cytokines found in {sample_name}, skipping...")
        continue
    
    # Find macrophage layers
    mac_layers, mac_obs = find_macrophage_layers(adata)
    logging.info(f"Found macrophage layers: {mac_layers}")
    logging.info(f"Found macrophage obs columns: {mac_obs}")
    
    # Create spatial plots for each present cytokine
    for gene in present:
        logging.info(f"Creating spatial plots for {gene}")
        
        # Plot from main expression matrix
        success = create_spatial_plot(adata, gene, sample_name=sample_name)
        
        # Plot from each macrophage layer if available
        for layer in mac_layers:
            success = create_spatial_plot(adata, gene, layer_name=layer, sample_name=sample_name)
    
    log_memory_usage()

logging.info("Immune secretion analysis complete!")
