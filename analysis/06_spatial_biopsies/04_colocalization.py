#!/usr/bin/env python3
# PROVENANCE: Ported 2026-06-05 from
#   /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/09_colocalization_analysis.py
# Generates the cell-type co-localization heatmap used as accepted-manuscript Fig 8I.
# Output filename (matched panel in bundle): mean_colocalization_heatmap_combined_macrophages.svg
# Input pickles:
#   /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeistNeilAnalysis/CITEgeist/analysis/data/biopsy_adatas.pkl
# Depends on utils.standardize_adata + figure_config.setup_figure_params from CITEgeist analysis dir.
import scanpy as sc
import os
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import gc
import logging
import json
import pickle  # Add this import for loading the pickle file
from sklearn.metrics.pairwise import cosine_similarity
from utils import standardize_adata # Assuming utils.py is in the Python path or same directory
from figure_config import setup_figure_params

# --- Configuration ---
# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Parse command line arguments
parser = argparse.ArgumentParser(description='Colocalization analysis for spatial transcriptomics')
parser.add_argument('--figures-only', action='store_true',
                    help='Skip computation; regenerate figures from saved colocalization matrices')
args = parser.parse_args()

# Apply standardized figure parameters (14pt minimum, Arial font)
setup_figure_params()

# Output directory for figures
output_dir = 'figures/colocalization_analysis'
os.makedirs(output_dir, exist_ok=True)
logging.info(f"Output directory set to: {output_dir}")

# Define standardized cell type names (based on cell_profiles.json)
standardized_celltypes = [
    'Cancer_Cells',
    'Immunosuppressive_Macrophages_CD163plus',
    'Inflammatory_Macrophages_CD14plusHLA_minusDRplus',
    'T_Cell–Interacting_Macrophages_CD11cplus',
    'CD4_T_Cells',
    'CD8_T_Cells',
    'B_Cells',
    'Endothelial_Cells',
    'Fibroblasts'
]

# Define combined cell types list (with merged macrophages)
combined_celltypes = [
    'Cancer_Cells',
    'Immunosuppressive_Macrophages_CD163plus',
    'Immune_Stimulating_Macrophages',  # This replaces the two separate macrophage types
    'CD4_T_Cells',
    'CD8_T_Cells',
    'B_Cells',
    'Endothelial_Cells',
    'Fibroblasts'
]

def calculate_colocalization_matrix(biopsy_adatas, cell_profiles, celltypes_list, output_name):
    """Calculate colocalization matrix for given cell types."""
    all_sample_coloc_matrices = []

    logging.info(f"Starting colocalization analysis for {output_name}...")
    for sample_name, adata in biopsy_adatas.items():
        # Skip specific sample types if needed
        if 'TNBC' in sample_name or 'NN165' in sample_name:
            logging.info(f"Skipping sample: {sample_name}")
            continue

        logging.info(f"Processing sample: {sample_name}")

        try:
            # Standardize adata.obs labels
            adata, name_mapping = standardize_adata(adata, cell_profiles)
            logging.info(f"Standardized AnnData for {sample_name}")

            # Filter celltypes_list to those present in this sample's adata.obs
            present_celltypes_in_sample = [ct for ct in celltypes_list if ct in adata.obs.columns]

            if len(present_celltypes_in_sample) < 2:
                logging.warning(f"Sample {sample_name} has fewer than 2 cell types with proportion data. Skipping colocalization calculation for this sample.")
                # Create a NaN matrix for this sample
                sample_coloc_df_reindexed = pd.DataFrame(np.nan, index=celltypes_list, columns=celltypes_list)
            else:
                # Extract the cell type proportions for present cell types
                proportions_df = adata.obs[present_celltypes_in_sample]

                # Calculate pairwise cosine similarity between cell types (columns)
                # cosine_similarity expects (n_samples, n_features).
                # To get similarity between columns (cell types), we transpose proportions_df
                # so that cell types become rows (n_celltypes) and spots become features (n_spots).
                similarity_values = cosine_similarity(proportions_df.T)

                # Create a DataFrame with proper labels
                sample_coloc_df = pd.DataFrame(similarity_values, index=present_celltypes_in_sample, columns=present_celltypes_in_sample)

                # Reindex to ensure all celltypes_list are present, filling with NaN
                # if a cell type was missing in this sample (already handled by present_celltypes_in_sample for calculation)
                # This ensures all matrices can be averaged correctly.
                sample_coloc_df_reindexed = sample_coloc_df.reindex(index=celltypes_list, columns=celltypes_list)

            all_sample_coloc_matrices.append(sample_coloc_df_reindexed)
            logging.info(f"Calculated colocalization matrix for {sample_name}")

        except Exception as e:
            logging.error(f"Error processing sample {sample_name}: {e}")
            # Append a NaN matrix in case of error to avoid breaking the aggregation
            nan_matrix = pd.DataFrame(np.nan, index=celltypes_list, columns=celltypes_list)
            all_sample_coloc_matrices.append(nan_matrix)
        finally:
            # Clean up memory
            gc.collect()

    logging.info("Finished processing all samples.")

    # --- Aggregate Colocalization Scores ---
    if not all_sample_coloc_matrices:
        logging.warning("No colocalization matrices were generated. Cannot create a mean heatmap.")
        mean_coloc_matrix = pd.DataFrame(np.nan, index=celltypes_list, columns=celltypes_list) # Empty plot
    else:
        logging.info(f"Aggregating {len(all_sample_coloc_matrices)} colocalization matrices...")
        
        # Check if there are any valid matrices to average
        valid_matrices = [m for m in all_sample_coloc_matrices if not m.isnull().all().all()]

        if not valid_matrices:
            logging.warning("All generated colocalization matrices are empty or NaN. Cannot compute mean.")
            mean_coloc_matrix = pd.DataFrame(np.nan, index=celltypes_list, columns=celltypes_list)
        else:
            # Use np.nanmean to average, ignoring NaNs. This is robust if some samples
            # didn't have all cell types or had errors.
            stacked_matrices = np.array([df.to_numpy() for df in valid_matrices])
            mean_coloc_values = np.nanmean(stacked_matrices, axis=0)
            mean_coloc_matrix = pd.DataFrame(mean_coloc_values, index=celltypes_list, columns=celltypes_list)
            logging.info("Successfully aggregated colocalization scores.")

    return mean_coloc_matrix

def plot_colocalization_heatmap(mean_coloc_matrix, title, filename):
    """Generate and save colocalization heatmap."""
    logging.info(f"Generating {title}...")
    plt.figure(figsize=(12, 10)) # Adjusted for potentially many cell types
    try:
        # Replace underscores with spaces in axis labels for better readability
        display_matrix = mean_coloc_matrix.copy()
        display_matrix.index = display_matrix.index.str.replace('_', ' ')
        display_matrix.columns = display_matrix.columns.str.replace('_', ' ')

        # Create a mask for the upper triangle to avoid redundant display (keep diagonal)
        mask = np.triu(np.ones_like(display_matrix, dtype=bool), k=1)

        ax = sns.heatmap(display_matrix, mask=mask, annot=True, cmap="viridis", fmt=".2f", vmin=0, vmax=1,
                    linewidths=0, linecolor='none', cbar_kws={'label': 'Mean Cosine Similarity'})
        # Remove edge colors and disable antialiasing to eliminate grid lines
        for patch in ax.collections:
            patch.set_edgecolor('face')  # Match edge to face color
            patch.set_antialiased(False)
        # Turn off axis grid lines
        ax.grid(False)
        plt.title(title, fontsize=16)
        plt.xticks(rotation=45, ha='right', fontsize=14)
        plt.yticks(rotation=0, fontsize=14)
        plt.tight_layout() # Adjust layout to prevent labels overlapping

        output_plot_path = os.path.join(output_dir, filename)
        # Use rasterized=True to reduce anti-aliasing artifacts at cell boundaries
        plt.savefig(output_plot_path, dpi=300, bbox_inches='tight', metadata={'Software': None})
        plt.savefig(output_plot_path.replace('.png', '.svg'), bbox_inches='tight')
        logging.info(f"Colocalization heatmap saved to {output_plot_path}")
    except Exception as e:
        logging.error(f"Error generating or saving heatmap: {e}")
    finally:
        plt.close() # Close the plot figure to free memory

# ============================================================================
# FIGURES-ONLY MODE: Skip computation, regenerate figures from saved matrices
# ============================================================================
if args.figures_only:
    logging.info("FIGURES-ONLY mode: loading saved colocalization matrices")

    coloc_pickle_path = 'data/colocalization_matrices.pkl'
    if not os.path.exists(coloc_pickle_path):
        logging.error(f"Cannot run figures-only mode: {coloc_pickle_path} not found")
        logging.error("Run full analysis first to generate intermediate results")
        sys.exit(1)

    with open(coloc_pickle_path, 'rb') as f:
        saved_matrices = pickle.load(f)

    original_coloc_matrix = saved_matrices['original']
    mean_coloc_matrix_combined = saved_matrices['combined']
    logging.info("Successfully loaded saved colocalization matrices")

    # Generate original heatmap
    plot_colocalization_heatmap(
        original_coloc_matrix,
        'Mean Cell Type Colocalization Across Biopsy Samples',
        'mean_colocalization_heatmap_original.png'
    )

    # Generate combined heatmap
    plot_colocalization_heatmap(
        mean_coloc_matrix_combined,
        'Mean Cell Type Colocalization Across Biopsy Samples',
        'mean_colocalization_heatmap_combined_macrophages.png'
    )

    logging.info("Figures-only mode complete!")
    sys.exit(0)

# ============================================================================
# FULL ANALYSIS MODE: Load data and run complete analysis
# ============================================================================

# Load cell profiles from JSON file
cell_profiles_path = 'cell_profiles.json'
try:
    with open(cell_profiles_path, 'r') as f:
        cell_profiles = json.load(f)
    logging.info(f"Successfully loaded cell profiles from {cell_profiles_path}")
except FileNotFoundError:
    logging.error(f"Error: {cell_profiles_path} not found. Please ensure the file exists.")
    sys.exit(1)
except json.JSONDecodeError:
    logging.error(f"Error: Could not decode JSON from {cell_profiles_path}.")
    sys.exit(1)

# Load biopsy data from pickle file
biopsy_adatas_path = 'data/biopsy_adatas.pkl'
try:
    with open(biopsy_adatas_path, 'rb') as f:
        biopsy_adatas = pickle.load(f)
    logging.info(f"Successfully loaded {len(biopsy_adatas)} biopsy samples from {biopsy_adatas_path}")
except FileNotFoundError:
    logging.error(f"Error: {biopsy_adatas_path} not found. Please ensure the file exists.")
    sys.exit(1)
except Exception as e:
    logging.error(f"Error loading pickle file: {e}")
    sys.exit(1)

# --- Original Analysis ---
# Calculate original colocalization matrix
original_coloc_matrix = calculate_colocalization_matrix(
    biopsy_adatas, cell_profiles, standardized_celltypes, "original cell types"
)

# Generate original heatmap
plot_colocalization_heatmap(
    original_coloc_matrix, 
    'Mean Cell Type Colocalization Across Biopsy Samples', 
    'mean_colocalization_heatmap_original.png'
)

# --- Combined Macrophage Analysis ---
logging.info("Starting combined macrophage analysis...")

# Create modified colocalization analysis with combined macrophages
all_sample_coloc_matrices_combined = []

for sample_name, adata in biopsy_adatas.items():
    # Skip specific sample types if needed
    if 'TNBC' in sample_name or 'NN165' in sample_name:
        logging.info(f"Skipping sample: {sample_name}")
        continue

    logging.info(f"Processing sample with combined macrophages: {sample_name}")

    try:
        # Standardize adata.obs labels
        adata, name_mapping = standardize_adata(adata, cell_profiles)
        
        # Create a modified observation dataframe with combined macrophages
        adata_obs_combined = adata.obs.copy()
        
        # Combine the two macrophage types into "Immune Stimulating Macrophages"
        inflammatory_col = 'Inflammatory_Macrophages_CD14plusHLA_minusDRplus'
        cd11c_col = 'T_Cell–Interacting_Macrophages_CD11cplus'
        combined_col = 'Immune_Stimulating_Macrophages'
        
        if inflammatory_col in adata_obs_combined.columns and cd11c_col in adata_obs_combined.columns:
            # Combine the proportions by taking the maximum of the two (since they represent different cell types)
            adata_obs_combined[combined_col] = adata_obs_combined[[inflammatory_col, cd11c_col]].sum(axis=1)
            
            # Remove the original columns from analysis
            adata_obs_combined = adata_obs_combined.drop(columns=[inflammatory_col, cd11c_col])
            logging.info(f"Combined {inflammatory_col} and {cd11c_col} into {combined_col}")
        elif inflammatory_col in adata_obs_combined.columns:
            # Only inflammatory present
            adata_obs_combined[combined_col] = adata_obs_combined[inflammatory_col]
            adata_obs_combined = adata_obs_combined.drop(columns=[inflammatory_col])
        elif cd11c_col in adata_obs_combined.columns:
            # Only CD11c+ present
            adata_obs_combined[combined_col] = adata_obs_combined[cd11c_col]
            adata_obs_combined = adata_obs_combined.drop(columns=[cd11c_col])
        else:
            # Neither present, add zero column
            adata_obs_combined[combined_col] = 0.0
            
        # Define combined cell types list
        combined_celltypes = [
            'Cancer_Cells',
            'Immunosuppressive_Macrophages_CD163plus',
            'Immune_Stimulating_Macrophages',  # This replaces the two separate macrophage types
            'CD4_T_Cells',
            'CD8_T_Cells',
            'B_Cells',
            'Endothelial_Cells',
            'Fibroblasts'
        ]
        
        # Filter combined_celltypes to those present in this sample's modified obs
        present_celltypes_in_sample = [ct for ct in combined_celltypes if ct in adata_obs_combined.columns]

        if len(present_celltypes_in_sample) < 2:
            logging.warning(f"Sample {sample_name} has fewer than 2 combined cell types with proportion data. Skipping colocalization calculation for this sample.")
            # Create a NaN matrix for this sample
            sample_coloc_df_reindexed = pd.DataFrame(np.nan, index=combined_celltypes, columns=combined_celltypes)
        else:
            # Extract the cell type proportions for present cell types
            proportions_df = adata_obs_combined[present_celltypes_in_sample]

            # Calculate pairwise cosine similarity between cell types (columns)
            similarity_values = cosine_similarity(proportions_df.T)

            # Create a DataFrame with proper labels
            sample_coloc_df = pd.DataFrame(similarity_values, index=present_celltypes_in_sample, columns=present_celltypes_in_sample)

            # Reindex to ensure all combined_celltypes are present, filling with NaN
            sample_coloc_df_reindexed = sample_coloc_df.reindex(index=combined_celltypes, columns=combined_celltypes)

        all_sample_coloc_matrices_combined.append(sample_coloc_df_reindexed)
        logging.info(f"Calculated combined colocalization matrix for {sample_name}")

    except Exception as e:
        logging.error(f"Error processing sample {sample_name} with combined macrophages: {e}")
        # Append a NaN matrix in case of error to avoid breaking the aggregation
        nan_matrix = pd.DataFrame(np.nan, index=combined_celltypes, columns=combined_celltypes)
        all_sample_coloc_matrices_combined.append(nan_matrix)
    finally:
        # Clean up memory
        gc.collect()

# --- Aggregate Combined Colocalization Scores ---
if not all_sample_coloc_matrices_combined:
    logging.warning("No combined colocalization matrices were generated. Cannot create a mean heatmap.")
    mean_coloc_matrix_combined = pd.DataFrame(np.nan, index=combined_celltypes, columns=combined_celltypes)
else:
    logging.info(f"Aggregating {len(all_sample_coloc_matrices_combined)} combined colocalization matrices...")
    
    # Check if there are any valid matrices to average
    valid_matrices_combined = [m for m in all_sample_coloc_matrices_combined if not m.isnull().all().all()]

    if not valid_matrices_combined:
        logging.warning("All generated combined colocalization matrices are empty or NaN. Cannot compute mean.")
        mean_coloc_matrix_combined = pd.DataFrame(np.nan, index=combined_celltypes, columns=combined_celltypes)
    else:
        # Use np.nanmean to average, ignoring NaNs
        stacked_matrices_combined = np.array([df.to_numpy() for df in valid_matrices_combined])
        mean_coloc_values_combined = np.nanmean(stacked_matrices_combined, axis=0)
        mean_coloc_matrix_combined = pd.DataFrame(mean_coloc_values_combined, index=combined_celltypes, columns=combined_celltypes)
        logging.info("Successfully aggregated combined colocalization scores.")

# Save intermediate results for figures-only mode
coloc_pickle_path = 'data/colocalization_matrices.pkl'
saved_matrices = {
    'original': original_coloc_matrix,
    'combined': mean_coloc_matrix_combined
}
with open(coloc_pickle_path, 'wb') as f:
    pickle.dump(saved_matrices, f)
logging.info(f"Saved colocalization matrices to {coloc_pickle_path}")

# Generate combined heatmap
plot_colocalization_heatmap(
    mean_coloc_matrix_combined, 
    'Mean Cell Type Colocalization Across Biopsy Samples', 
    'mean_colocalization_heatmap_combined_macrophages.png'
)

logging.info("Script finished.")
