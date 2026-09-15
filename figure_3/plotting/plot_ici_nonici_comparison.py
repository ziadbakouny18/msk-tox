"""Render fig3b/fig3c ICI vs Non-ICI AE comparison figures from aggregate CSVs.

Port of ``plotting/plot_ici_nonici_comparison.py`` (plus the bar/label half of
``plotting/plot_ici_nonici_barplot.py``) with all computation removed. C-index
values, bootstrap CIs, KM cumulative-incidence step tables, multivariate
logrank p-values, and per-panel y-max are read from the source-data CSVs; the
script only draws.

Panel sets:
    A (fig3b): 5 grade0_base AEs
    B (fig3c): pan-cancer all-AE, grade 0
    C (fig3c): pan-cancer all-AE, grade 3
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter, NullLocator

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
try:
    from aesthetics import set_nature_style
except ImportError:  # imported as a module instead of run as a script
    import sys

    if _SCRIPT_DIR not in sys.path:
        sys.path.insert(0, _SCRIPT_DIR)
    from aesthetics import set_nature_style

_ICI_COLOR = "#3C5488"
_NON_ICI_COLOR = "#E64B35"
_PANEL_WIDTH_INCHES = 0.9
_ROW_HEIGHT_UNIT_INCHES = 0.72
_HORIZONTAL_GAP_INCHES = 0.4
_VERTICAL_GAP_INCHES = 0.18
_OUTER_PANEL_GAP_INCHES = _HORIZONTAL_GAP_INCHES / 2
_LEFT_MARGIN_INCHES = 0.72
_RIGHT_MARGIN_INCHES = 0.18
_TOP_MARGIN_INCHES = 0.52
_BOTTOM_MARGIN_INCHES = 0.42
_FONT_SIZE_PT = 7
_RISK_ORDER = ["High", "Medium", "Low"]
_RISK_COLORS = {
    "High": "#c0392b",
    "Medium": "#e67e22",
    "Low": "#2980b9",
}
_OK_STATUS = "ok"
_SUBSETS = ("ICI", "Non-ICI")
_PANEL_SET_FILENAMES = {
    "A": "fig3b_ici_nonici_comparison",
    "B": "fig3c_ici_nonici_comparison",
    "C": "fig3c_ici_nonici_comparison_grade3",
}
_DEFAULT_SOURCE_DATA = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "source_data"))
_DEFAULT_OUT_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "figures"))


@dataclass
class _Panel:
    ae: str
    order: int
    ymax: float
    cindex: dict[str, pd.Series]
    logrank: dict[str, pd.Series]
    curves: dict[str, dict[str, pd.DataFrame]]


def _pretty_panel_label(ae: str) -> str:
    return "All AEs" if str(ae).strip().lower() == "all" else str(ae).replace("_", " ").title()


def _load_sources(source_data: str, panel_set: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    def read(name: str) -> pd.DataFrame:
        path = os.path.join(source_data, name)
        df = pd.read_csv(path)
        if "panel_set" not in df.columns:
            raise ValueError(f"panel_set column missing in {path}")
        return df.loc[df["panel_set"] == panel_set].reset_index(drop=True)

    meta = read("ici_nonici_panel_meta.csv")
    cindex = read("ici_nonici_cindex.csv")
    logrank = read("ici_nonici_logrank.csv")
    curves = read("ici_nonici_incidence_curves.csv")
    return meta, cindex, logrank, curves


def _build_panels(meta: pd.DataFrame, cindex: pd.DataFrame, logrank: pd.DataFrame, curves: pd.DataFrame) -> list[_Panel]:
    panels: list[_Panel] = []
    for _, row in meta.sort_values("panel_order").iterrows():
        ae = str(row["ae"])
        c_rows = cindex.loc[cindex["ae"] == ae].set_index("subset")
        l_rows = logrank.loc[logrank["ae"] == ae].set_index("subset")
        sub = curves.loc[curves["ae"] == ae]
        curve_map: dict[str, dict[str, pd.DataFrame]] = {}
        for subset in _SUBSETS:
            s2 = sub.loc[sub["subset"] == subset]
            curve_map[subset] = {
                str(risk): g.sort_values("time_days").reset_index(drop=True)
                for risk, g in s2.groupby("risk_tertile")
            }
        panels.append(
            _Panel(
                ae=ae,
                order=int(row["panel_order"]),
                ymax=float(row["ymax"]),
                cindex={s: c_rows.loc[s] for s in _SUBSETS if s in c_rows.index},
                logrank={s: l_rows.loc[s] for s in _SUBSETS if s in l_rows.index},
                curves=curve_map,
            )
        )
    return panels


def _max_time_days(curves: pd.DataFrame) -> int:
    times = pd.to_numeric(curves.get("time_days"), errors="coerce").dropna()
    return int(times.max()) if not times.empty else 365


def _draw_insufficient_panel(ax: Axes, n_samples: int, n_events: int) -> None:
    ax.text(0.5, 0.5, f"Insufficient data\nN={n_samples}, events={n_events}", ha="center", va="center", fontsize=_FONT_SIZE_PT, transform=ax.transAxes)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])


def _bar_value(row: pd.Series | None) -> float | None:
    if row is None or row["status"] != _OK_STATUS or pd.isna(row["c_index"]):
        return None
    return float(row["c_index"])


def _plot_cindex_bars(ax: Axes, panel: _Panel, show_ylabel: bool) -> None:
    width = 0.36
    values: dict[str, float | None] = {}
    for subset, pos in (("ICI", -width / 2), ("Non-ICI", width / 2)):
        row = panel.cindex.get(subset)
        value = _bar_value(row)
        values[subset] = value
        if value is None:
            continue
        yerr = np.zeros((2, 1), dtype=float)
        if pd.notna(row["ci_low"]) and pd.notna(row["ci_high"]):
            yerr = np.array([[max(0.0, value - float(row["ci_low"]))], [max(0.0, float(row["ci_high"]) - value)]], dtype=float)
        ax.bar([pos], [value], width=width, color=_ICI_COLOR if subset == "ICI" else _NON_ICI_COLOR, yerr=yerr, capsize=2, error_kw={"elinewidth": 0.7, "capthick": 0.7, "ecolor": "#2f2f2f"})
    if values["ICI"] is None and values["Non-ICI"] is None:
        ax.text(0.5, 0.02, "NA", ha="center", va="bottom", fontsize=_FONT_SIZE_PT, transform=ax.transAxes)
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0.5, 0.9)
    ax.set_xticks([])
    ax.set_yticks(np.arange(0.5, 0.91, 0.1))
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    ax.tick_params(axis="y", which="both", direction="out", left=True, right=False, labelsize=_FONT_SIZE_PT)
    ax.yaxis.set_ticks_position("left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel("C-index" if show_ylabel else "")


def _month_tick_positions(max_time_days: int) -> tuple[np.ndarray, list[str]]:
    if max_time_days <= 0:
        return np.array([0.0]), ["0"]

    max_months = max(1, int(round((12 * max_time_days) / 365.0)))
    step_months = 3 if max_months >= 6 else 1
    month_values = np.arange(0, max_months + 1, step_months, dtype=int)
    if month_values[-1] != max_months:
        month_values = np.append(month_values, max_months)

    tick_positions = month_values.astype(float) * (max_time_days / max_months)
    tick_positions[-1] = float(max_time_days)
    return tick_positions, [str(int(month)) for month in month_values]


def _km_subset_ok(panel: _Panel, subset: str) -> bool:
    row = panel.logrank.get(subset)
    return row is not None and row["status"] == _OK_STATUS


def _subset_counts(panel: _Panel, subset: str) -> tuple[int, int]:
    row = panel.cindex.get(subset)
    n_samples = int(row["n_samples"]) if row is not None and pd.notna(row["n_samples"]) else 0
    n_events = int(row["n_events"]) if row is not None and pd.notna(row["n_events"]) else 0
    return n_samples, n_events


def _plot_quantile_cumulative_incidence(
    ax: Axes,
    panel: _Panel,
    subset: str,
    panel_title: str,
    max_time_days: int,
    y_upper: float,
    show_ylabel: bool,
    show_xlabel: bool,
) -> None:
    if not _km_subset_ok(panel, subset):
        n_samples, n_events = _subset_counts(panel, subset)
        _draw_insufficient_panel(ax, n_samples, n_events)
        ax.set_ylim(0, y_upper)
        ax.set_title(panel_title, fontsize=_FONT_SIZE_PT, fontweight="bold")
        return

    for risk in _RISK_ORDER:
        steps = panel.curves[subset].get(risk)
        if steps is None or steps.empty:
            continue
        steps = steps.dropna(subset=["time_days", "cum_incidence"])
        if steps.empty:
            continue
        color = _RISK_COLORS[risk]
        t = steps["time_days"].to_numpy(dtype=float)
        y = steps["cum_incidence"].to_numpy(dtype=float)
        if not steps["ci_low"].isna().all() and not steps["ci_high"].isna().all():
            ax.fill_between(t, steps["ci_low"].to_numpy(dtype=float), steps["ci_high"].to_numpy(dtype=float), step="post", color=color, alpha=0.12, linewidth=0)
        ax.step(t, y, where="post", color=color, linewidth=1.0)

    logrank_row = panel.logrank[subset]
    p_val = float(logrank_row["p_value"])
    p_text = f"p = {p_val:.3f}" if p_val >= 0.001 else "p < 0.001"
    n_samples, n_events = _subset_counts(panel, subset)
    ax.text(0.03, 0.97, f"{p_text}\nLoT = {n_samples}\nAE = {n_events}", transform=ax.transAxes, ha="left", va="top", fontsize=_FONT_SIZE_PT)

    ax.set_title(panel_title, fontsize=_FONT_SIZE_PT, fontweight="bold")
    ax.set_ylabel("Incidence" if show_ylabel else "")
    if show_xlabel:
        ax.set_xlabel("Months")
    else:
        ax.set_xlabel("")

    ax.set_xlim(0, max_time_days)
    ax.set_ylim(0, y_upper)
    xticks, xticklabels = _month_tick_positions(max_time_days)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(
        axis="x",
        which="major",
        direction="out" if show_xlabel else "in",
        bottom=True,
        top=False,
        labelbottom=show_xlabel,
        labelsize=_FONT_SIZE_PT,
    )
    ax.tick_params(axis="x", which="minor", bottom=False, top=False)

    ystep = 0.05 if y_upper < 0.199 else 0.1
    ax.set_yticks(np.arange(0, y_upper + ystep * 0.5, ystep))
    ax.tick_params(axis="y", which="both", direction="out", left=True, right=False, labelsize=_FONT_SIZE_PT)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f" if ystep < 0.1 else "%.1f"))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _plot_three_row_figure(panels: list[_Panel], out_prefix: str, max_time_days: int) -> None:
    set_nature_style()
    plt.rcParams.update({"font.size": _FONT_SIZE_PT, "axes.labelsize": _FONT_SIZE_PT, "axes.titlesize": _FONT_SIZE_PT, "xtick.labelsize": _FONT_SIZE_PT, "ytick.labelsize": _FONT_SIZE_PT, "legend.fontsize": _FONT_SIZE_PT, "legend.title_fontsize": _FONT_SIZE_PT})

    n_ae = len(panels)
    row_height_ratios = [1.0, 1.25, 1.25]
    panel_heights = [ratio * _ROW_HEIGHT_UNIT_INCHES for ratio in row_height_ratios]
    total_panel_width = (n_ae * _PANEL_WIDTH_INCHES) + max(0, n_ae - 1) * _HORIZONTAL_GAP_INCHES
    inner_width = total_panel_width + 2 * _OUTER_PANEL_GAP_INCHES
    fig_width = _LEFT_MARGIN_INCHES + inner_width + _RIGHT_MARGIN_INCHES
    total_panel_height = sum(panel_heights) + 2 * _VERTICAL_GAP_INCHES
    fig_height = _TOP_MARGIN_INCHES + total_panel_height + _BOTTOM_MARGIN_INCHES
    fig = plt.figure(figsize=(fig_width, fig_height))
    average_panel_height = sum(panel_heights) / len(panel_heights)
    gs = GridSpec(
        3,
        ncols=n_ae,
        figure=fig,
        height_ratios=row_height_ratios,
        wspace=_HORIZONTAL_GAP_INCHES / _PANEL_WIDTH_INCHES,
        hspace=_VERTICAL_GAP_INCHES / average_panel_height,
    )
    fig.subplots_adjust(
        left=(_LEFT_MARGIN_INCHES + _OUTER_PANEL_GAP_INCHES) / fig_width,
        right=1 - ((_RIGHT_MARGIN_INCHES + _OUTER_PANEL_GAP_INCHES) / fig_width),
        top=1 - (_TOP_MARGIN_INCHES / fig_height),
        bottom=_BOTTOM_MARGIN_INCHES / fig_height,
    )

    ref_ax: Axes | None = None
    first_ici_ax: Axes | None = None
    first_non_ax: Axes | None = None
    for i, panel in enumerate(panels):
        ax_bar = fig.add_subplot(gs[0, i]) if i == 0 else fig.add_subplot(gs[0, i], sharey=ref_ax)
        if i == 0:
            ref_ax = ax_bar
        _plot_cindex_bars(ax_bar, panel, show_ylabel=(i == 0))

    for i, panel in enumerate(panels):
        ae_label = _pretty_panel_label(panel.ae)
        ax_ici = fig.add_subplot(gs[1, i])
        _plot_quantile_cumulative_incidence(ax_ici, panel, "ICI", ae_label, max_time_days, panel.ymax, show_ylabel=(i == 0), show_xlabel=False)
        if i == 0:
            first_ici_ax = ax_ici

        ax_non = fig.add_subplot(gs[2, i], sharex=ax_ici)
        _plot_quantile_cumulative_incidence(ax_non, panel, "Non-ICI", "", max_time_days, panel.ymax, show_ylabel=(i == 0), show_xlabel=True)
        for shared_ax in (ax_ici, ax_non):
            shared_ax.xaxis.set_minor_locator(NullLocator())
            shared_ax.tick_params(axis="x", which="minor", bottom=False, top=False)
        if i == 0:
            first_non_ax = ax_non

    if first_ici_ax is not None:
        ici_box = first_ici_ax.get_position()
        fig.text(
            (_LEFT_MARGIN_INCHES * 0.32) / fig_width,
            ici_box.y0 + (ici_box.height / 2),
            "ICI",
            color=_ICI_COLOR,
            rotation=90,
            ha="center",
            va="center",
            fontsize=_FONT_SIZE_PT,
            fontweight="bold",
        )
    if first_non_ax is not None:
        non_box = first_non_ax.get_position()
        fig.text(
            (_LEFT_MARGIN_INCHES * 0.32) / fig_width,
            non_box.y0 + (non_box.height / 2),
            "Non-ICI",
            color=_NON_ICI_COLOR,
            rotation=90,
            ha="center",
            va="center",
            fontsize=_FONT_SIZE_PT,
            fontweight="bold",
        )

    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=_ICI_COLOR, edgecolor="black", linewidth=0.5, label="ICI"),
        Rectangle((0, 0), 1, 1, facecolor=_NON_ICI_COLOR, edgecolor="black", linewidth=0.5, label="Non-ICI"),
    ] + [Line2D([0], [0], color=_RISK_COLORS[risk], linewidth=1.0, label=risk) for risk in _RISK_ORDER]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        frameon=False,
        ncol=len(legend_handles),
        borderaxespad=0,
        title=None,
        fontsize=_FONT_SIZE_PT,
        handlelength=1.4,
        columnspacing=1.0,
    )
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(f"{out_prefix}.png", dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(f"{out_prefix}.pdf", dpi=600, facecolor="white", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render fig3b/fig3c ICI vs Non-ICI comparison figures from aggregate source-data CSVs.")
    parser.add_argument("--panel_set", required=True, choices=["A", "B", "C"], help="A = fig3b (5 grade0_base AEs), B = fig3c (all-AE grade0), C = fig3c (all-AE grade3).")
    parser.add_argument("--source_data", type=str, default=_DEFAULT_SOURCE_DATA, help="Directory containing the source-data CSVs (default: <figure_3>/source_data).")
    parser.add_argument("--out_dir", type=str, default=_DEFAULT_OUT_DIR, help="Directory for rendered PNG/PDF (default: <figure_3>/figures).")
    args = parser.parse_args()

    meta, cindex, logrank, curves = _load_sources(args.source_data, args.panel_set)
    panels = _build_panels(meta, cindex, logrank, curves)
    if not panels:
        raise RuntimeError(f"No panels found in source data for panel_set={args.panel_set}.")

    max_time_days = _max_time_days(curves)
    out_prefix = os.path.join(args.out_dir, _PANEL_SET_FILENAMES[args.panel_set])
    _plot_three_row_figure(panels, out_prefix, max_time_days)
    print(f"Saved figure: {out_prefix}.png")
    print(f"Saved figure: {out_prefix}.pdf")


if __name__ == "__main__":
    main()
