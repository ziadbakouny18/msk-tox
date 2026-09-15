"""Nature Publishing Group aesthetic definitions for MSK-Tox figures.

Vendored from OncoAE ``rsf/src/aesthetics.py`` so the figure-replication
package is self-contained. Import ``set_nature_style`` before plotting.
"""

from __future__ import annotations

import matplotlib as mpl
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


# ---------------------------------------------------------------------------
# Palettes
# ---------------------------------------------------------------------------

# Categorical palette (10 colours) following NPG conventions.
NATURE_COLORS: list[str] = [
    "#DC0000",  # bright red
    "#3C5488",  # blue
    "#E64B35",  # red-orange
    "#00A087",  # teal
    "#4DBBD5",  # cyan
    "#F39B7F",  # salmon
    "#8491B4",  # lavender
    "#91D1C2",  # mint
    "#7E6148",  # brown
    "#B09C85",  # tan
]

# Ordered palette for multi-model colour strips.
MODEL_TYPE_COLORS: list[str] = ["#3C5488", "#E64B35", "#00A087", "#F39B7F"]


# ---------------------------------------------------------------------------
# rcParams helper
# ---------------------------------------------------------------------------

def set_nature_style() -> None:
    """Configure matplotlib for Nature publication style."""
    sns.set_style("white")
    sns.set_context("paper", font_scale=1.0)
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica"]
    mpl.rcParams["font.size"] = 6
    mpl.rcParams["axes.labelsize"] = 7
    mpl.rcParams["axes.titlesize"] = 8
    mpl.rcParams["xtick.labelsize"] = 6
    mpl.rcParams["ytick.labelsize"] = 6
    mpl.rcParams["xtick.direction"] = "out"
    mpl.rcParams["ytick.direction"] = "out"
    mpl.rcParams["legend.fontsize"] = 5
    mpl.rcParams["legend.title_fontsize"] = 6
    mpl.rcParams["axes.linewidth"] = 0.5
    mpl.rcParams["xtick.major.width"] = 0.5
    mpl.rcParams["ytick.major.width"] = 0.5
    mpl.rcParams["xtick.major.size"] = 2
    mpl.rcParams["ytick.major.size"] = 2
    mpl.rcParams["lines.linewidth"] = 0.75
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["savefig.dpi"] = 600
    mpl.rcParams["savefig.transparent"] = False
    mpl.rcParams["savefig.facecolor"] = "white"
    mpl.rcParams["savefig.edgecolor"] = "none"
    mpl.rcParams["figure.facecolor"] = "white"
    mpl.rcParams["axes.facecolor"] = "white"
    mpl.rcParams["axes.edgecolor"] = "black"
    mpl.rcParams["axes.labelcolor"] = "black"
    mpl.rcParams["xtick.color"] = "black"
    mpl.rcParams["ytick.color"] = "black"
    mpl.rcParams["text.color"] = "black"


# ---------------------------------------------------------------------------
# Palette helpers
# ---------------------------------------------------------------------------

def get_nature_palette(n_colors: int, palette_type: str = "categorical") -> list:
    """Return Nature Publishing Group color palettes as RGB tuples.

    Args:
        n_colors: Number of colours requested.
        palette_type: ``"categorical"`` — NPG categorical palette (cycles if
            n_colors > 10). ``"model"`` — ordered colours for model-type
            strips. ``"sequential"`` — viridis sequential palette.
    """
    if palette_type == "model":
        return [mpl.colors.to_rgb(c) for c in MODEL_TYPE_COLORS[:n_colors]]
    elif palette_type == "sequential":
        return sns.color_palette("viridis", n_colors)
    else:
        colors = NATURE_COLORS[-n_colors:] if n_colors <= len(NATURE_COLORS) else NATURE_COLORS
        return [mpl.colors.to_rgb(c) for c in colors]


# ---------------------------------------------------------------------------
# Colormap
# ---------------------------------------------------------------------------

def get_diverging_cmap() -> LinearSegmentedColormap:
    """Colorblind-friendly blue-white-red diverging colormap."""
    colors = [
        "#2166AC", "#4393C3", "#92C5DE",
        "#D1E5F0", "#F7F7F7", "#FDDBC7",
        "#F4A582", "#D6604D", "#B2182B",
    ]
    return LinearSegmentedColormap.from_list("colorblind_diverging", colors)
