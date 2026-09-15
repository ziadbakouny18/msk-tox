"""Render supp3e AE-level c-index swarm plot across model configurations from aggregate CSVs.

Port of ``plotting/plot_cindex_swarm.py`` with all computation removed.
Per-AE, per-model c-index values, the paired-bootstrap p-value vs base, and
the FDR-corrected q-value are read from the source-data CSV; the script only
draws. Point significance (black outline) is taken from the
``p_value_vs_base_fdr`` column. The deterministic swarm jitter /
collision-offset layout from the original is preserved.
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
try:
    from aesthetics import NATURE_COLORS, set_nature_style
except ImportError:  # imported as a module instead of run as a script
    import sys

    if _SCRIPT_DIR not in sys.path:
        sys.path.insert(0, _SCRIPT_DIR)
    from aesthetics import NATURE_COLORS, set_nature_style

MODEL_ORDER = [
    ("base", ["base"], "base"),
    ("+labs", ["base", "labs"], "base+labs"),
    ("+tumor\ngenomics", ["base", "genomics"], "base+genomics"),
    ("+tumor\ngenomics\n+labs", ["base", "labs", "genomics"], "base+labs+genomics"),
]
_SIGNIFICANCE_THRESHOLD = 0.05
_CSV_FILE = "cindex_swarm.csv"
_OUTPUT_PREFIX = "supp3e_cindex_swarm"
_DEFAULT_SOURCE_DATA = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "source_data"))
_DEFAULT_OUT_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "figures"))


def _load_source(source_data: str) -> pd.DataFrame:
    path = os.path.join(source_data, _CSV_FILE)
    df = pd.read_csv(path)
    for column in ("c_index", "delta_vs_base", "p_value_vs_base", "p_value_vs_base_fdr", "n_samples", "n_events", "paired_n"):
        df[column] = pd.to_numeric(df[column], errors="coerce")
    return df.reset_index(drop=True)


def _swarm_positions(y_values: list[float], x_center: float, x_step: float = 0.045) -> list[float]:
    if not y_values:
        return []

    y_min = min(y_values)
    y_max = max(y_values)
    y_threshold = max(0.004, 0.06 * max(y_max - y_min, 0.02))
    placed: list[tuple[float, float]] = []
    positions = [x_center] * len(y_values)
    candidate_levels = [0]
    for level in range(1, 20):
        candidate_levels.extend([-level, level])

    for idx, y_value in sorted(enumerate(y_values), key=lambda item: item[1]):
        x_value = x_center
        for level in candidate_levels:
            candidate = x_center + level * x_step
            collision = any(
                abs(y_value - other_y) < y_threshold and abs(candidate - other_x) < x_step * 0.98
                for other_x, other_y in placed
            )
            if not collision:
                x_value = candidate
                break
        positions[idx] = x_value
        placed.append((x_value, y_value))
    return positions


def _make_palette(ae_names: list[str]) -> dict[str, str]:
    return {ae: NATURE_COLORS[idx % len(NATURE_COLORS)] for idx, ae in enumerate(ae_names)}


def _is_significant(p_value: float | None, alpha: float = _SIGNIFICANCE_THRESHOLD) -> bool:
    return p_value is not None and not np.isnan(p_value) and p_value < alpha


def _plot(summary: pd.DataFrame, output_prefix: str) -> None:
    set_nature_style()

    x_labels = [label for label, _, _ in MODEL_ORDER]
    x_pos = {label: idx for idx, label in enumerate(x_labels)}
    ae_names = sorted(str(ae) for ae in summary.loc[summary["c_index"].notna(), "ae"].unique())
    palette = _make_palette(ae_names)

    plot_df = summary[summary["c_index"].notna()].copy()
    plot_df["x_center"] = plot_df["model_label"].map(x_pos).astype(float)
    plot_df["x_plot"] = plot_df["x_center"]
    plot_df["is_significant"] = plot_df["p_value_vs_base_fdr"].apply(_is_significant)
    for model_label in x_labels:
        mask = plot_df["model_label"] == model_label
        y_values = plot_df.loc[mask, "c_index"].tolist()
        plot_df.loc[mask, "x_plot"] = _swarm_positions(y_values, x_center=float(x_pos[model_label]))

    fig, ax = plt.subplots(figsize=(3, 3))

    for ae, df_ae in plot_df.sort_values("x_center").groupby("ae"):
        ae_name = str(ae)
        ax.plot(
            df_ae["x_plot"],
            df_ae["c_index"],
            color=palette[ae_name],
            linewidth=0.8,
            alpha=0.55,
            zorder=2,
        )

    for ae, df_ae in plot_df.groupby("ae"):
        ae_name = str(ae)
        edgecolors = ["black" if is_significant else "none" for is_significant in df_ae["is_significant"]]
        linewidths = [1.5 if is_significant else 0.0 for is_significant in df_ae["is_significant"]]
        ax.scatter(
            df_ae["x_plot"],
            df_ae["c_index"],
            s=20,
            color=palette[ae_name],
            edgecolors=edgecolors,
            linewidths=linewidths,
            zorder=3,
            label=ae_name,
        )

    c_values = plot_df["c_index"].to_numpy(dtype=float)
    y_min = 0.5
    y_max = min(1.0, float(np.nanmax(c_values)) + 0.10)
    ax.set_ylim(y_min, y_max)

    ax.tick_params(axis="both", which="both", direction="out", width=1.0, length=3)
    ax.set_xticks(list(x_pos.values()))
    ax.set_xticklabels(x_labels)
    ax.set_xlim(-0.5, len(x_labels) - 0.5)
    ax.set_ylabel("c-index")
    ax.set_xlabel("")
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color=palette[ae],
            markerfacecolor=palette[ae],
            markeredgecolor="black",
            markeredgewidth=0.35,
            markersize=4,
            linewidth=0.8,
            label=ae.replace("_", " "),
        )
        for ae in ae_names
    ]

    significance_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#BDBDBD",
            markeredgecolor="black",
            markeredgewidth=0.6,
            markersize=5,
            linewidth=0,
            label="q<0.05 vs base",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#BDBDBD",
            markeredgecolor="none",
            markeredgewidth=0.0,
            markersize=5,
            linewidth=0,
            label="q>=0.05 vs base",
        ),
    ]

    legend_significance = ax.legend(
        handles=significance_handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
        ncol=1,
    )
    ax.add_artist(legend_significance)
    ax.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.02, 0.42), frameon=False, ncol=1)

    fig.subplots_adjust(left=0.12, right=0.78, bottom=0.18, top=0.96)

    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(f"{output_prefix}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{output_prefix}.pdf", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render supp3e AE-level c-index swarm plot across model configurations from aggregate source-data CSVs."
    )
    parser.add_argument("--source_data", type=str, default=_DEFAULT_SOURCE_DATA, help="Directory containing the source-data CSVs (default: <figure_3>/source_data).")
    parser.add_argument("--out_dir", type=str, default=_DEFAULT_OUT_DIR, help="Directory for rendered PNG/PDF (default: <figure_3>/figures).")
    args = parser.parse_args()

    summary = _load_source(args.source_data)
    if summary["c_index"].notna().sum() == 0:
        raise ValueError("No valid c-index values were found in cindex_swarm.csv.")

    output_prefix = os.path.join(args.out_dir, _OUTPUT_PREFIX)
    _plot(summary, output_prefix)
    print(f"Saved figure: {output_prefix}.png")
    print(f"Saved figure: {output_prefix}.pdf")


if __name__ == "__main__":
    main()
