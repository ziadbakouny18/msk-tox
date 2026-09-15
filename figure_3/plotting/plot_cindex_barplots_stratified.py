"""Render supp3a/supp3b c-index barplots stratified by LOT or treatment modality from aggregate CSVs.

Port of ``plotting/plot_cindex_barplots_stratified.py`` with all computation
removed. Per-AE, per-stratum c-index values, bootstrap CIs, and the
unpaired-bootstrap p-value vs the baseline stratum are read from the
source-data CSVs; the script only draws.

Panels:
    supp3a (--stratify_by lot):        5 grade0_base AEs x line-of-therapy strata
    supp3b (--stratify_by treatment):  5 grade0_base AEs x treatment-modality strata
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
try:
    from aesthetics import set_nature_style
except ImportError:  # imported as a module instead of run as a script
    import sys

    if _SCRIPT_DIR not in sys.path:
        sys.path.insert(0, _SCRIPT_DIR)
    from aesthetics import set_nature_style

STRATIFIED_LOT_GROUP_ORDER = ["Line 1", "Line 2", "Line 3", "Line 4+"]
TREATMENT_GROUP_ORDER = [
    "PD-(L)1",
    "CTLA-4",
    "Biologic",
    "Hormonal",
    "Chemotherapy",
    "Targeted",
    "Other",
]
# AE display order matches the figure; the CSVs list AEs alphabetically.
AE_ORDER = ["adrenal_insufficiency", "hypothyroidism", "liver_toxicity", "pneumonitis", "colitis"]
MIN_TREATMENT_PATIENTS = 20
HALF_LETTER_WIDTH_INCHES = 4.25
SIGNIFICANCE_THRESHOLD = 0.05
SIGNIFICANCE_BRACKET_OFFSET = 0.02
SIGNIFICANCE_BRACKET_HEIGHT = 0.009
SIGNIFICANCE_LABEL_GAP = 0.004
SIGNIFICANCE_BRACKET_GAP = 0.04

_CSV_FILES = {
    "lot": "cindex_stratified_lot.csv",
    "treatment": "cindex_stratified_treatment.csv",
}
_OUTPUT_PREFIXES = {
    "lot": "supp3a_cindex_stratified_lot",
    "treatment": "supp3b_cindex_stratified_treatment",
}
_DEFAULT_SOURCE_DATA = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "source_data"))
_DEFAULT_OUT_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "figures"))


def _load_source(source_data: str, stratify_by: str) -> pd.DataFrame:
    path = os.path.join(source_data, _CSV_FILES[stratify_by])
    df = pd.read_csv(path).rename(
        columns={
            "stratum": "group_label",
            "ci_low": "bootstrap_ci_low",
            "ci_high": "bootstrap_ci_high",
            "p_value": "p_value_vs_first",
        }
    )
    if "included_in_figure" in df.columns:
        flag = df["included_in_figure"].astype(str).str.strip().str.lower()
        df = df[flag.isin(["true", "1", "yes"])]
    return df.reset_index(drop=True)


def _build_palette(group_order: list[str], stratify_by: str) -> dict[str, tuple[float, float, float, float]]:
    if stratify_by == "lot":
        colors = ["#c6dbef", "#6baed6", "#2171b5", "#08519c"]
    elif stratify_by == "treatment":
        cmap = plt.get_cmap("Set2")
        colors = [cmap(i) for i in range(len(group_order))]
    else:
        raise ValueError(f"Unknown stratification mode: {stratify_by}")
    return {group: colors[i] for i, group in enumerate(group_order)}


def _fmt_p_sig(p: float | None) -> str:
    if p is None or np.isnan(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def _draw_significance_bracket(ax, x_left: float, x_right: float, y: float, p_value: float) -> float:
    bracket_top = y + SIGNIFICANCE_BRACKET_HEIGHT
    ax.plot(
        [x_left, x_left, x_right, x_right],
        [y, bracket_top, bracket_top, y],
        color="black",
        linewidth=0.7,
        clip_on=False,
    )
    ax.text(
        (x_left + x_right) / 2.0,
        bracket_top + SIGNIFICANCE_LABEL_GAP,
        _fmt_p_sig(p_value),
        ha="center",
        va="bottom",
        fontsize=5,
        color="black",
        clip_on=False,
    )
    return bracket_top + SIGNIFICANCE_LABEL_GAP


def _pretty_ae_label(ae: str) -> str:
    return str(ae).replace("_", " ").title()


def _plot_summary(summary_df: pd.DataFrame, *, ae_names: list[str], stratify_by: str, output_prefix: str) -> None:
    set_nature_style()
    if stratify_by == "treatment":
        n_samples = pd.to_numeric(summary_df["n_samples"], errors="coerce")
        eligible_groups = set(summary_df.loc[n_samples >= MIN_TREATMENT_PATIENTS, "group_label"].astype(str))
        summary_df = summary_df[summary_df["group_label"].astype(str).isin(eligible_groups)]
    group_order = list(summary_df["group_label"].dropna().astype(str).unique())
    group_order = [group for group in (STRATIFIED_LOT_GROUP_ORDER if stratify_by == "lot" else TREATMENT_GROUP_ORDER) if group in group_order]
    if not group_order:
        raise RuntimeError(f"No groups available to plot for mode: {stratify_by}")
    palette = _build_palette(group_order, stratify_by)
    n_groups = len(group_order)
    width = min(0.8 / max(1, n_groups), 0.22)
    x = np.arange(len(ae_names))
    fig, ax = plt.subplots(1, 1, figsize=(HALF_LETTER_WIDTH_INCHES, 3.0))
    indexed = summary_df.set_index(["ae", "group_label"]).sort_index()
    ci_tops = np.full(len(ae_names), 0.5, dtype=float)
    bar_centers: dict[str, np.ndarray] = {}
    bar_tops: dict[str, np.ndarray] = {}
    significant_comparisons: dict[int, list[tuple[str, float]]] = {index: [] for index in range(len(ae_names))}
    baseline_groups: dict[str, str] = {}
    for ae in ae_names:
        baseline = summary_df.loc[summary_df["ae"] == ae, "baseline_group"].dropna()
        if not baseline.empty:
            baseline_groups[ae] = str(baseline.iloc[0])

    for i, group in enumerate(group_order):
        offset = (i - (n_groups - 1) / 2.0) * width
        x_group = x + offset
        y_vals: list[float] = []
        yerr_low: list[float] = []
        yerr_high: list[float] = []
        p_vals: list[float | None] = []
        dropped_by_patient_threshold: list[bool] = []
        for ae in ae_names:
            key = (ae, group)
            if key in indexed.index:
                row = indexed.loc[key]
                n_samples = pd.to_numeric(row["n_samples"], errors="coerce")
                drop_bar = stratify_by == "treatment" and (pd.isna(n_samples) or n_samples < MIN_TREATMENT_PATIENTS)
                c_val = row["c_index"]
                ci_low = row["bootstrap_ci_low"]
                ci_high = row["bootstrap_ci_high"]
                p_val = row["p_value_vs_first"]
                value = np.nan if drop_bar or pd.isna(c_val) else float(c_val)
                y_vals.append(value)
                if not np.isnan(value) and pd.notna(ci_low) and pd.notna(ci_high):
                    yerr_low.append(max(0.0, value - float(ci_low)))
                    yerr_high.append(max(0.0, float(ci_high) - value))
                else:
                    yerr_low.append(0.0)
                    yerr_high.append(0.0)
                p_vals.append(None if pd.isna(p_val) else float(p_val))
                dropped_by_patient_threshold.append(drop_bar)
            else:
                y_vals.append(np.nan)
                yerr_low.append(0.0)
                yerr_high.append(0.0)
                p_vals.append(None)
                dropped_by_patient_threshold.append(stratify_by == "treatment")

        bars = ax.bar(x_group, np.asarray(y_vals, dtype=float), width=width, color=palette[group], label=group)
        bar_centers[group] = np.asarray([bar.get_x() + bar.get_width() / 2.0 for bar in bars], dtype=float)
        bar_tops[group] = np.asarray(
            [value + error_high if not np.isnan(value) else np.nan for value, error_high in zip(y_vals, yerr_high)],
            dtype=float,
        )
        for j, bar in enumerate(bars):
            y_val = y_vals[j]
            if np.isnan(y_val):
                if not dropped_by_patient_threshold[j]:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2.0,
                        0.515,
                        "NA",
                        ha="center",
                        va="bottom",
                        fontsize=4.5,
                        color="#808080",
                    )
                continue
            ci_tops[j] = max(ci_tops[j], y_val + yerr_high[j])
            p_val = p_vals[j]
            errorbar_color = "black" if p_val is not None and p_val < SIGNIFICANCE_THRESHOLD else "#808080"
            ax.errorbar(
                bar.get_x() + bar.get_width() / 2.0,
                y_val,
                yerr=np.array([[yerr_low[j]], [yerr_high[j]]]),
                fmt="none",
                capsize=2,
                elinewidth=0.7,
                capthick=0.7,
                ecolor=errorbar_color,
            )
            if p_val is not None and p_val < SIGNIFICANCE_THRESHOLD:
                significant_comparisons[j].append((group, p_val))

    bracket_text_tops: list[float] = []
    for ae_index, comparisons in significant_comparisons.items():
        baseline_group = baseline_groups.get(ae_names[ae_index])
        if baseline_group is None or baseline_group not in bar_centers:
            continue
        baseline_x = bar_centers[baseline_group][ae_index]
        baseline_top = bar_tops[baseline_group][ae_index]
        if np.isnan(baseline_top):
            continue
        next_bracket_y = -np.inf
        for comparison_group, p_value in comparisons:
            if comparison_group not in bar_centers:
                continue
            comparison_top = bar_tops[comparison_group][ae_index]
            if np.isnan(comparison_top):
                continue
            comparison_x = bar_centers[comparison_group][ae_index]
            bracket_y = max(baseline_top, comparison_top) + SIGNIFICANCE_BRACKET_OFFSET
            bracket_y = max(bracket_y, next_bracket_y)
            x_left, x_right = sorted((baseline_x, comparison_x))
            text_top = _draw_significance_bracket(ax, x_left, x_right, bracket_y, p_value)
            bracket_text_tops.append(text_top)
            next_bracket_y = text_top + SIGNIFICANCE_BRACKET_GAP

    ax.set_xticks(x)
    ax.set_xticklabels([_pretty_ae_label(ae) for ae in ae_names], rotation=25, ha="right")
    max_annotation_y = max(
        max((ci_tops[ae_index] + SIGNIFICANCE_BRACKET_OFFSET for ae_index in range(len(ae_names)) if ci_tops[ae_index] > 0.5), default=1.0),
        max(bracket_text_tops, default=1.0),
    )
    ax.set_ylim(0.4, max(1.03, max_annotation_y + 0.015))
    ax.set_ylabel("c-index")
    ax.tick_params(axis="y", direction="out")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=len(group_order), frameon=False, fontsize=4.5, columnspacing=0.8, handletextpad=0.3)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.80, bottom=0.28)
    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(f"{output_prefix}.png", dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{output_prefix}.pdf", dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Render supp3a/supp3b c-index barplots stratified by LOT or treatment modality from aggregate source-data CSVs.")
    parser.add_argument("--stratify_by", type=str, choices=["lot", "treatment"], required=True, help="Stratification mode: lot (supp3a) or treatment (supp3b).")
    parser.add_argument("--source_data", type=str, default=_DEFAULT_SOURCE_DATA, help="Directory containing the source-data CSVs (default: <figure_3>/source_data).")
    parser.add_argument("--out_dir", type=str, default=_DEFAULT_OUT_DIR, help="Directory for rendered PNG/PDF (default: <figure_3>/figures).")
    args = parser.parse_args()

    summary_df = _load_source(args.source_data, args.stratify_by)
    ae_names = [ae for ae in AE_ORDER if ae in set(summary_df["ae"].astype(str))]
    if not ae_names:
        raise RuntimeError(f"No AEs available to plot for mode: {args.stratify_by}.")

    output_prefix = os.path.join(args.out_dir, _OUTPUT_PREFIXES[args.stratify_by])
    _plot_summary(summary_df, ae_names=ae_names, stratify_by=args.stratify_by, output_prefix=output_prefix)
    print(f"Saved figure: {output_prefix}.png")
    print(f"Saved figure: {output_prefix}.pdf")


if __name__ == "__main__":
    main()
