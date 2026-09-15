"""Render supp3c SHAP beeswarm plots for Renal Cell Carcinoma (RCC) per AE from a long-format source CSV.

Port of the SHAP summary (beeswarm) rendering chain in
``rsf/src/shap_analysis.py`` (``render_shap_plots``) with all computation
removed. Per-dot SHAP values, feature values, importance ranks, binary
flags, and display labels are read from the source-data CSV; the script
only draws. The CSV rows were independently permuted within each
(ae, feature) group for de-identification, so per-row dot y-jitter may
differ slightly from the original renders.
"""

from __future__ import annotations

import argparse
import os
from typing import Any, cast

import numpy as np
import pandas as pd
import shap
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PathCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
try:
    from aesthetics import set_nature_style
except ImportError:  # imported as a module instead of run as a script
    import sys

    if _SCRIPT_DIR not in sys.path:
        sys.path.insert(0, _SCRIPT_DIR)
    from aesthetics import set_nature_style

# Vendored color constants from ``rsf/src/shap_analysis.py``.
_CONTINUOUS_LOW_COLOR = "#4C6C86"
_CONTINUOUS_MID_COLOR = "#C8C1C6"
_CONTINUOUS_HIGH_COLOR = "#A0556B"
_BINARY_ZERO_COLOR = "#5D7D7A"
_BINARY_ONE_COLOR = "#A3604F"
_MISSING_FEATURE_COLOR = "#B5B5B5"
_CONTINUOUS_CMAP = LinearSegmentedColormap.from_list(
    "oncoae_shap_continuous",
    [_CONTINUOUS_LOW_COLOR, _CONTINUOUS_MID_COLOR, _CONTINUOUS_HIGH_COLOR],
)

_MAX_DISPLAY = 15
_CSV_FILE = "shap_beeswarm_rcc_shuffled_deidentified.csv"
_OUTPUT_PREFIX = "supp3c_shap_beeswarm"
_AE_CHOICES = ("liver_toxicity", "hypothyroidism")
_DEFAULT_SOURCE_DATA = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "source_data"))
_DEFAULT_OUT_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "figures"))


def _pretty_title_token(value: str) -> str:
    return str(value).replace("_", " ").strip().title()


def _load_source(source_data: str, ae: str) -> pd.DataFrame:
    path = os.path.join(source_data, _CSV_FILE)
    df = pd.read_csv(path)
    df = df[df["ae"] == ae].copy()
    if df.empty:
        raise ValueError(f"No rows for ae={ae!r} in {path}.")
    return df.reset_index(drop=True)


def _build_plot_inputs(
    df: pd.DataFrame,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, dict[str, bool]]:
    """Pivot the long CSV to (n_dots x n_features) shap values + feature values.

    Returns the shap matrix, the display-label DataFrame passed to
    ``shap.summary_plot``, the raw-name DataFrame, and the per-feature
    binary flag. Columns are ordered by ascending ``feature_rank``
    (rank 0 = top row); ``shap.summary_plot`` re-sorts rows by mean |SHAP|
    internally, so the drawn layout matches ``feature_rank``.
    """
    order = (
        df.groupby("feature")["feature_rank"].min().sort_values(kind="stable").index.tolist()
    )
    counts = df.groupby("feature")["shap_value"].count()
    n = int(counts[order[0]])
    if not (counts[order] == n).all():
        raise ValueError("Unequal dot counts across features; cannot pivot to a matrix.")

    shap_cols: list[np.ndarray] = []
    value_cols: list[np.ndarray] = []
    labels: list[str] = []
    binary_flags: dict[str, bool] = {}
    for feature in order:
        group = df[df["feature"] == feature]
        shap_cols.append(group["shap_value"].to_numpy(dtype=float))
        value_cols.append(group["feature_value"].to_numpy(dtype=float))
        labels.append(str(group["feature_label"].iloc[0]))
        binary_flags[feature] = bool(group["is_binary"].iloc[0])

    shap_values = np.column_stack(shap_cols)
    X_display = pd.DataFrame(dict(zip(labels, value_cols)))
    X_raw = pd.DataFrame(dict(zip(order, value_cols)))
    return shap_values, X_display, X_raw, binary_flags


# ---------------------------------------------------------------------------
# Styling helpers vendored from ``rsf/src/shap_analysis.py``.
# ---------------------------------------------------------------------------


def _display_feature_order(shap_values: np.ndarray, max_display: int) -> np.ndarray:
    """Return feature indices in the same bottom-to-top order used by SHAP."""
    mean_abs = np.abs(shap_values).mean(axis=0)
    top_feature_idx = np.argsort(mean_abs)[-int(max_display):]
    return top_feature_idx


def _row_collections(summary_ax: Axes, n_rows: int) -> dict[int, list[PathCollection]]:
    """Group non-empty SHAP scatter collections by their plotted feature row."""
    grouped: dict[int, list[PathCollection]] = {row: [] for row in range(n_rows)}

    for collection in summary_ax.collections:
        if not isinstance(collection, PathCollection):
            continue
        offsets = np.asarray(collection.get_offsets())
        if offsets.size == 0:
            continue

        y_median = float(np.nanmedian(offsets[:, 1]))
        row_idx = int(round(y_median))
        if 0 <= row_idx < n_rows:
            grouped[row_idx].append(collection)

    return grouped


def _set_collection_constant_color(
    collection: PathCollection,
    color: str,
    alpha: float,
) -> None:
    offsets = np.asarray(collection.get_offsets())
    rgba = np.tile(to_rgba(color, alpha=alpha), (len(offsets), 1))
    collection_any = cast(Any, collection)
    collection.set_array(None)
    collection_any.set_facecolors(rgba)
    collection.set_edgecolor("none")


def _recolor_shap_summary_rows(
    summary_ax: Axes,
    X_plot: pd.DataFrame,
    shap_values: np.ndarray,
    max_display: int,
    binary_flags: dict[str, bool],
) -> tuple[bool, bool]:
    """Apply separate muted color treatments to binary and continuous rows.

    Port of ``_recolor_shap_summary_rows`` from ``rsf/src/shap_analysis.py``.
    The binary/continuous classification is read from the CSV
    ``is_binary`` column instead of the feature name.
    """
    feature_order = _display_feature_order(shap_values, max_display=max_display)
    grouped_collections = _row_collections(summary_ax, n_rows=len(feature_order))

    has_binary = False
    has_continuous = False

    for row_idx, feature_idx in enumerate(feature_order):
        col_name = X_plot.columns[feature_idx]
        is_binary = binary_flags[col_name]
        row_has_data = False

        for collection in grouped_collections.get(row_idx, []):
            row_has_data = True
            feature_array = collection.get_array()

            if is_binary:
                if feature_array is None:
                    _set_collection_constant_color(collection, _MISSING_FEATURE_COLOR, alpha=0.50)
                else:
                    values = np.asarray(feature_array, dtype=float)
                    rgba = np.tile(
                        to_rgba(_BINARY_ZERO_COLOR, alpha=0.58),
                        (len(values), 1),
                    )
                    rgba[values >= 0.5] = to_rgba(_BINARY_ONE_COLOR, alpha=0.58)
                    collection_any = cast(Any, collection)
                    collection.set_array(None)
                    collection_any.set_facecolors(rgba)
                    collection.set_edgecolor("none")
                has_binary = True
            else:
                if feature_array is None:
                    _set_collection_constant_color(collection, _MISSING_FEATURE_COLOR, alpha=0.50)
                else:
                    values = np.asarray(feature_array, dtype=float)
                    finite_values = values[np.isfinite(values)]
                    if finite_values.size:
                        vmin, vmax = np.nanpercentile(finite_values, [5, 95])
                        if not np.isfinite(vmin) or not np.isfinite(vmax) or np.isclose(vmin, vmax):
                            vmin = float(np.nanmin(finite_values))
                            vmax = float(np.nanmax(finite_values))
                            if np.isclose(vmin, vmax):
                                vmax = vmin + 1.0
                        collection.set_clim(vmin, vmax)
                    collection.set_cmap(_CONTINUOUS_CMAP)
                    collection.set_edgecolor("none")
                has_continuous = True

        if not row_has_data and is_binary:
            has_binary = True
        elif not row_has_data and not is_binary:
            has_continuous = True

    return has_binary, has_continuous


def _add_shap_summary_legends(
    fig: Figure,
    summary_ax: Axes,
    has_binary: bool,
    has_continuous: bool,
) -> None:
    """Replace the default single SHAP legend with vertically stacked
    binary + continuous legends.

    Layout (in figure coordinates, sharing the same left edge):

        Binary value          <- legend title (horizontal)
          ● 0
          ● 1

        Continuous value      <- colorbar title (horizontal, same style)
          ┃ High
          ┃  …   (short vertical colorbar)
          ┃ Low
    """
    legend_left = 0.76             # shared left edge for both blocks
    colorbar_ax = fig.axes[-1] if len(fig.axes) > 1 else None

    marker_size = 5.0              # legend dot diameter, in points
    # -- Binary value legend (top block) ------------------------------------
    if has_binary:
        legend_handles = [
            Line2D(
                [0], [0], marker="o", linestyle="None", markersize=marker_size,
                markerfacecolor=_BINARY_ZERO_COLOR, markeredgecolor="none",
                alpha=0.72, label="False",
            ),
            Line2D(
                [0], [0], marker="o", linestyle="None", markersize=marker_size,
                markerfacecolor=_BINARY_ONE_COLOR, markeredgecolor="none",
                alpha=0.72, label="True",
            ),
        ]
        binary_legend = fig.legend(
            handles=legend_handles,
            title="Binary value",
            loc="upper left",
            bbox_to_anchor=(legend_left, 0.88),
            frameon=False,
            fontsize=4.5,
            title_fontsize=5.5,
            handletextpad=0.35,
            labelspacing=0.35,
            borderaxespad=0.0,
        )
        binary_legend._legend_box.align = "left"
        fig.add_artist(binary_legend)

    # -- Continuous value colorbar (bottom block, stacked below binary) -----
    if colorbar_ax is not None:
        colorbar_ax.clear()
        if has_continuous:
            # Short vertical colorbar placed directly beneath the binary block
            # (or higher up when there is no binary legend above it).
            bar_bottom = 0.50 if has_binary else 0.56
            bar_height = 0.22 if has_binary else 0.28
            # Match the bar thickness to the binary legend dot diameter.
            # A colorbar enforces a default ~20:1 aspect that would otherwise
            # shrink the width we set, so we release it and re-assert geometry
            # AFTER the colorbar is built.
            fig_width_pt = fig.get_size_inches()[0] * 72.0
            bar_width = (marker_size * 1.25) / fig_width_pt
            bar_rect = (legend_left + 0.006, bar_bottom, bar_width, bar_height)
            colorbar_ax.set_position(bar_rect)

            scalar_mappable = ScalarMappable(norm=Normalize(vmin=0.0, vmax=1.0), cmap=_CONTINUOUS_CMAP)
            scalar_mappable.set_array([])
            colorbar = fig.colorbar(scalar_mappable, cax=colorbar_ax, orientation="vertical")
            colorbar.set_ticks([0.0, 1.0])
            colorbar.set_ticklabels(["Low", "High"])
            cast(Any, colorbar.outline).set_visible(False)
            # Release the colorbar aspect lock so the width above is honored.
            colorbar_ax.set_aspect("auto")
            colorbar_ax.set_box_aspect(None)
            colorbar_ax.set_position(bar_rect)
            # Horizontal title above the colorbar, matching the "Binary value" style.
            colorbar_ax.set_title(
                "Continuous value",
                fontsize=5.5,
                loc="left",
                pad=4,
            )
            colorbar_ax.tick_params(axis="y", which="both", labelsize=4.5, length=0, pad=1.5)
        else:
            colorbar_ax.set_visible(False)


def _restyle_shap_summary_figure(fig: Figure, summary_ax: Axes) -> None:
    """Tighten SHAP beeswarm text hierarchy and marker styling."""
    fig.patch.set_facecolor("white")
    summary_ax.set_facecolor("white")
    summary_ax.set_title(summary_ax.get_title(), fontsize=4.5, pad=3)
    summary_ax.set_xlabel(summary_ax.get_xlabel(), fontsize=5.5, labelpad=3)
    summary_ax.set_ylabel(summary_ax.get_ylabel(), fontsize=5.5, labelpad=4)
    summary_ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        top=False,
        right=False,
        labelsize=4.5,
        pad=1.5,
    )

    for collection in summary_ax.collections:
        if isinstance(collection, PathCollection):
            collection.set_sizes(np.full(collection.get_sizes().shape, 3.5))
            collection.set_alpha(0.8)
            collection.set_edgecolor("none")
            collection.set_rasterized(True)

    for axis in fig.axes:
        axis.set_facecolor("white")
        axis.tick_params(axis="both", which="both", direction="out", labelsize=4.5, pad=1.5)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


def _plot(
    shap_values: np.ndarray,
    X_display: pd.DataFrame,
    X_raw: pd.DataFrame,
    binary_flags: dict[str, bool],
    title_str: str,
    output_prefix: str,
) -> None:
    set_nature_style()
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica"],
            "font.size": 4.5,
            "axes.labelsize": 5.5,
            "axes.titlesize": 4.5,
            "xtick.labelsize": 4.5,
            "ytick.labelsize": 4.5,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 2.0,
            "ytick.major.size": 2.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    n_display = min(len(X_raw.columns), _MAX_DISPLAY)
    summary_height = max(4.4, min(5.8, 0.21 * n_display + 1.6))
    summary_width = 6.8

    # shap.summary_plot uses the current figure; do NOT call plt.figure()
    # before it or the manually-created figure leaks.
    shap.summary_plot(
        shap_values,
        X_display,
        max_display=_MAX_DISPLAY,
        show=False,
        plot_size=(summary_width, summary_height),
    )
    summary_fig = plt.gcf()
    summary_ax = plt.gca()
    summary_ax.set_title(f"{title_str}", fontsize=4.5, pad=3)
    has_binary, has_continuous = _recolor_shap_summary_rows(
        summary_ax=summary_ax,
        X_plot=X_raw,
        shap_values=shap_values,
        max_display=_MAX_DISPLAY,
        binary_flags=binary_flags,
    )
    _restyle_shap_summary_figure(summary_fig, summary_ax)
    # Set directional x-axis label AFTER restyle (restyle re-applies get_xlabel()).
    summary_ax.set_xlabel(
        "SHAP value (impact on predicted risk)", fontsize=5.5, labelpad=3
    )
    # Directional end-cues clarify the cumulative-hazard sign convention:
    # positive SHAP → higher predicted risk / earlier AE onset.
    _cue_style = dict(
        transform=summary_ax.transAxes,
        fontsize=4.5,
        color="#555555",
        va="top",
        clip_on=False,
    )
    summary_ax.text(0.0, -0.07, "← lower risk / later onset", ha="left", **_cue_style)
    summary_ax.text(1.0, -0.07, "higher risk / earlier onset →", ha="right", **_cue_style)
    summary_fig.subplots_adjust(left=0.34, right=0.82, top=0.90, bottom=0.16)
    _add_shap_summary_legends(
        fig=summary_fig,
        summary_ax=summary_ax,
        has_binary=has_binary,
        has_continuous=has_continuous,
    )
    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    summary_path = f"{output_prefix}.pdf"
    plt.savefig(summary_path, bbox_inches="tight", facecolor="white")
    summary_png_path = f"{output_prefix}.png"
    plt.savefig(summary_png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(summary_fig)
    print(f"  SHAP summary saved: {summary_path}")
    print(f"  SHAP summary saved: {summary_png_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render supp3c SHAP beeswarm plots for RCC per AE from the long-format source-data CSV."
    )
    parser.add_argument("--ae", type=str, choices=_AE_CHOICES, required=True, help="Adverse event to render.")
    parser.add_argument("--source_data", type=str, default=_DEFAULT_SOURCE_DATA, help="Directory containing the source-data CSV (default: <figure_3>/source_data).")
    parser.add_argument("--out_dir", type=str, default=_DEFAULT_OUT_DIR, help="Directory for rendered PNG/PDF (default: <figure_3>/figures).")
    args = parser.parse_args()

    df = _load_source(args.source_data, args.ae)
    shap_values, X_display, X_raw, binary_flags = _build_plot_inputs(df)
    cancer_type = str(df["cancer_type"].iloc[0])
    title_str = f"{_pretty_title_token(args.ae)} in {_pretty_title_token(cancer_type)}"
    output_prefix = os.path.join(args.out_dir, f"{_OUTPUT_PREFIX}_{args.ae}_rcc")
    _plot(shap_values, X_display, X_raw, binary_flags, title_str, output_prefix)


if __name__ == "__main__":
    main()
