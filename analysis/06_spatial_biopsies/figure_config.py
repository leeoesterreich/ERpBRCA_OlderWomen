"""Centralized figure configuration for spatial biopsy analysis scripts.

This module standardizes figure generation across the project with:
- Minimum 14pt font size for all text elements
- Arial font family (with fallback)
- Consistent DPI and figure settings
- Dual-format (SVG + PNG) save helper
"""

import os
import matplotlib.pyplot as plt
from matplotlib import font_manager
import scanpy as sc
from pathlib import Path
from typing import Optional, Union
import logging

# Font configuration — with fallback chain
ARIAL_FONT_PATH = '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Arial.ttf'
if not os.path.exists(ARIAL_FONT_PATH):
    # Fallback to system Liberation Sans (metrically equivalent to Arial)
    ARIAL_FONT_PATH = '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf'
if not os.path.exists(ARIAL_FONT_PATH):
    ARIAL_FONT_PATH = None  # matplotlib will use default
MIN_FONT_SIZE = 14

# Standard rcParams for all figures
FIGURE_RCPARAMS = {
    # Font settings
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'legend.title_fontsize': 14,

    # Figure settings
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',

    # Layout settings
    'figure.autolayout': False,  # Use explicit tight_layout instead
    'figure.constrained_layout.use': False,
}


def setup_figure_params(
    fontsize: int = MIN_FONT_SIZE,
    use_arial: bool = True,
    scanpy_params: bool = True
) -> None:
    """Apply standardized figure parameters.

    Args:
        fontsize: Base font size (minimum 14pt enforced)
        use_arial: Whether to load and use Arial font
        scanpy_params: Whether to also configure scanpy figure params
    """
    # Enforce minimum font size
    fontsize = max(fontsize, MIN_FONT_SIZE)

    # Load Arial font if available
    if use_arial:
        try:
            if Path(ARIAL_FONT_PATH).exists():
                font_manager.fontManager.addfont(ARIAL_FONT_PATH)
                plt.rcParams['font.family'] = 'Arial'
        except Exception as e:
            logging.warning(f"Could not load Arial font: {e}")

    # Apply base rcParams
    plt.rcParams.update(FIGURE_RCPARAMS)

    # Override with scaled font sizes
    plt.rcParams.update({
        'font.size': fontsize,
        'axes.labelsize': fontsize + 2,
        'axes.titlesize': fontsize + 4,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'legend.fontsize': fontsize,
        'legend.title_fontsize': fontsize,
    })

    # Configure scanpy if requested
    if scanpy_params:
        sc.set_figure_params(
            scanpy=True,
            fontsize=fontsize,
            dpi=100,
            dpi_save=300,
            frameon=False,
            figsize=(10, 10)
        )


def save_figure(
    fig_or_path: Union[plt.Figure, str, Path],
    output_path: Union[str, Path],
    dpi: int = 300,
    tight: bool = True,
    close: bool = True
) -> None:
    """Save figure in both PNG and SVG formats with tight layout.

    Args:
        fig_or_path: matplotlib Figure object or path string for plt.savefig
        output_path: Base output path (without extension, or with .png/.svg)
        dpi: DPI for PNG output
        tight: Whether to use bbox_inches='tight'
        close: Whether to close the figure after saving
    """
    output_path = Path(output_path)

    # Remove extension if present to create base path
    if output_path.suffix in ['.png', '.svg', '.pdf']:
        base_path = output_path.with_suffix('')
    else:
        base_path = output_path

    png_path = base_path.with_suffix('.png')
    svg_path = base_path.with_suffix('.svg')

    # Determine if we have a figure object or should use current figure
    if isinstance(fig_or_path, plt.Figure):
        fig = fig_or_path
    else:
        fig = plt.gcf()

    # Apply tight layout if requested
    if tight:
        try:
            fig.tight_layout()
        except Exception:
            pass  # Some figures don't support tight_layout

    # Save in both formats
    bbox = 'tight' if tight else None
    fig.savefig(png_path, dpi=dpi, bbox_inches=bbox)
    fig.savefig(svg_path, bbox_inches=bbox)

    if close:
        plt.close(fig)


def get_annotation_fontsize() -> int:
    """Return fontsize for heatmap annotations (minimum 14pt)."""
    return MIN_FONT_SIZE


def get_text_fontsize() -> int:
    """Return fontsize for general text annotations."""
    return MIN_FONT_SIZE


# Default figure sizes adjusted for larger fonts
FIGURE_SIZES = {
    'small': (8, 6),
    'medium': (10, 8),
    'large': (12, 10),
    'wide': (14, 8),
    'square': (10, 10),
    'heatmap': (12, 10),
    'spatial': (10, 10),
}
