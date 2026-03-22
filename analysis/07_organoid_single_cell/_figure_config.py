"""Shared figure configuration for NeilOrganoidSingleCell analysis.

Import this module early (after matplotlib.use('Agg')) to set:
- Arial font from custom .ttf
- Publication-quality font sizes
- save_fig() helper that outputs both PNG and SVG
"""
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Register Arial font
ARIAL_PATH = '/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/Arial.ttf'
if os.path.exists(ARIAL_PATH):
    fm.fontManager.addfont(ARIAL_PATH)
    matplotlib.rcParams['font.family'] = 'Arial'
else:
    matplotlib.rcParams['font.family'] = 'sans-serif'

# Global rcParams: Arial font, larger sizes
plt.rcParams.update({
    'font.family': matplotlib.rcParams['font.family'],
    'font.sans-serif': ['Arial'],
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 18,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
    'figure.constrained_layout.use': True,
})


def save_fig(path, **kwargs):
    """Save current figure as both PNG and SVG.

    Parameters
    ----------
    path : str
        Output path (should end in .png). SVG saved alongside.
    **kwargs
        Extra kwargs passed to plt.savefig (dpi, bbox_inches already defaulted).
    """
    kwargs.setdefault('dpi', 150)
    kwargs.setdefault('bbox_inches', 'tight')
    plt.savefig(path, **kwargs)
    svg_path = str(path).rsplit('.', 1)[0] + '.svg'
    plt.savefig(svg_path, format='svg', bbox_inches='tight')
