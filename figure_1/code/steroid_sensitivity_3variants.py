#!/usr/bin/env python3
"""
Steroid exposure sensitivity analysis — three metric definitions side by side.

Reads steroid_ratio_cache_LOT.pkl and renders the same cohort under three
definitions so they can be compared directly:

  V1  ratio, zeros excluded   log2((post+e)/(pre+e)), patients with any steroid
                              exposure in either window. This is the current
                              published Figure 1F-H definition.
  V2  ratio, zeros included   same value, but patients with no steroid in either
                              window are retained. They land at exactly 0.0
                              because (0+e)/(0+e) == 1.
  V3  post only               log2(post) among patients with post > 0. No pre
                              window, so no reference-window comparability
                              problem and no epsilon.

V3 discards patients with no post-window steroid rather than flooring them at
log2(e). That information is not thrown away: the proportion of patients with
any post-window steroid is reported per grade in the companion CSV, which is the
honest way to use the zeros (a proportion, not an imputed ratio).

Grade bins are FIXED at No AE / Grade 1-2 / Grade 3-4 for every toxicity and
every variant. The published notebook lets bin membership depend on per-panel
counts, which would change between variants and make them incomparable.

Outputs (results/sensitivity/):
  Steroid_Sensitivity_3variants.pdf / .png   3 variants x 3 toxicities
  Steroid_Sensitivity_3variants_stats.csv    per-bin stats for every variant
  Steroid_Sensitivity_by_raw_grade.csv       ungrouped grade 0-4 breakdown
"""

import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import kendalltau, mannwhitneyu, norm

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
ROOT = Path(__file__).parent.parent.resolve()
FIGURES = ROOT.parent.parent
DATA = FIGURES / "figures_data" / "figure 1" / "data"
CACHE_PATH = DATA / "steroid_ratio_cache_LOT.pkl"
OUT_DIR = ROOT / "results" / "sensitivity"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Must match generate_steroid_cache.py so V2 reproduces the cached log2_ratio.
EPSILON = 0.01
ICI_COL = "ici_at_ref"

TOXICITIES = [
    ("adrenal_insufficiency", "Adrenal Insufficiency"),
    ("colitis", "Colitis"),
    ("pneumonitis", "Pneumonitis"),
]

# Fixed bins, identical across variants and toxicities.
BINS = [("No AE", [0]), ("Grade 1-2", [1, 2]), ("Grade 3-4", [3, 4])]

GRADE_GRAY = "#95a5a6"


def bin_colors():
    cmap = sns.color_palette("rocket_r", as_cmap=True)
    return [GRADE_GRAY, cmap(0.35), cmap(0.75)]


def jonckheere_terpstra(groups):
    """Jonckheere-Terpstra trend test (same implementation as the notebook).

    WARNING: this counts strictly-greater pairs and uses the no-ties variance.
    Both assumptions fail when a large share of values are exactly equal, which
    is the case for V2 (46% of values sit at exactly 0.0 because a patient with
    no steroid in either window gets (0+e)/(0+e) == 1). Tied pairs contribute
    zero to the statistic while the null mean still assumes they split evenly,
    so z is deflated toward — or past — zero. Use trend_tau_b as the primary
    result; this is retained only for continuity with the published notebook.
    """
    k = len(groups)
    n_total = sum(len(g) for g in groups)
    JT = 0
    for i in range(k):
        for j in range(i + 1, k):
            gi, gj = np.asarray(groups[i]), np.asarray(groups[j])
            if len(gi) and len(gj):
                JT += int((gj[:, None] > gi[None, :]).sum())
    n_i = np.array([len(g) for g in groups])
    mean_JT = (n_total**2 - np.sum(n_i**2)) / 4
    var_JT = (n_total**2 * (2 * n_total + 3) - np.sum(n_i**2 * (2 * n_i + 3))) / 72
    if var_JT <= 0:
        return JT, np.nan, np.nan
    z = (JT - mean_JT) / np.sqrt(var_JT)
    log_p = np.log(2) + norm.logsf(abs(z))
    return JT, z, (np.exp(log_p) if log_p > -700 else 0.0)


def trend_tau(groups):
    """Tie-corrected monotonic trend: Kendall tau-b of bin index vs value.

    Directionally equivalent to Jonckheere-Terpstra but handles ties correctly
    in both the ordinal predictor and the outcome, and tau-b is an effect size
    (bounded -1..1, does not inflate with n). Returns (tau_b, p, pct_tied_at_0).
    """
    vals = np.concatenate([g for g in groups if len(g)]) if any(len(g) for g in groups) else np.array([])
    idx = np.concatenate([np.full(len(g), i) for i, g in enumerate(groups) if len(g)])
    if len(vals) < 3 or len(np.unique(idx)) < 2:
        return np.nan, np.nan, np.nan
    res = kendalltau(idx, vals, variant="b")
    pct_tied = 100.0 * float(np.mean(np.isclose(vals, 0.0)))
    return float(res.statistic), float(res.pvalue), pct_tied


def ptxt(p):
    if p is None or not np.isfinite(p):
        return "n/a"
    return "p<0.001" if p < 0.001 else f"p={p:.3f}"


# ---------------------------------------------------------------------------
# LOAD
# ---------------------------------------------------------------------------
with open(CACHE_PATH, "rb") as f:
    _raw = pickle.load(f)
if isinstance(_raw, dict) and _raw.get("_compat") == "numpy_columns_v1":
    _data = {}
    for _col in _raw["columns"]:
        _arr = _raw["arrays"][_col]
        _kind = _raw["kinds"][_col]
        if _kind == "datetime_ns":
            _data[_col] = pd.to_datetime(_arr, unit="ns")
        else:
            _data[_col] = _arr
    full = pd.DataFrame(_data, columns=_raw["columns"])
else:
    full = _raw

print(f"Cache: {CACHE_PATH}")
print(f"  {len(full):,} patient-toxicity rows")

ici = full[full[ICI_COL]].copy()
print(f"  ICI+ at reference date: {len(ici):,}")

# Recompute the ratio for every ICI+ row, including the all-zero patients that
# the cache stores as NaN. Reproduces log2_ratio exactly where it is defined.
ici["ratio_all"] = np.log2((ici["post_steroid"] + EPSILON) / (ici["pre_steroid"] + EPSILON))
defined = ici["log2_ratio"].notna()
assert np.allclose(ici.loc[defined, "ratio_all"], ici.loc[defined, "log2_ratio"]), (
    "recomputed ratio does not match cached log2_ratio -- check EPSILON"
)
print(f"  EPSILON check passed (eps={EPSILON})")

ici["log2_post"] = np.log2(ici["post_steroid"].where(ici["post_steroid"] > 0))


# ---------------------------------------------------------------------------
# VARIANTS
# ---------------------------------------------------------------------------
VARIANTS = [
    {
        "key": "V1_ratio_zeros_excluded",
        "row_label": "V1  ratio, zeros excluded\n(current Figure 1F-H)",
        "ylabel": "log2(post/pre)\nprednisone equiv. mg-days",
        "value": "ratio_all",
        "mask": lambda d: d["any_steroid"],
    },
    {
        "key": "V2_ratio_zeros_included",
        "row_label": "V2  ratio, zeros included\n(no-steroid patients kept at 0)",
        "ylabel": "log2(post/pre)\nprednisone equiv. mg-days",
        "value": "ratio_all",
        "mask": lambda d: pd.Series(True, index=d.index),
    },
    {
        "key": "V3_post_only",
        "row_label": "V3  post window only\n(no pre-window reference)",
        "ylabel": "log2(post)\nprednisone equiv. mg-days",
        "value": "log2_post",
        "mask": lambda d: d["post_steroid"] > 0,
    },
]


def collect(sub, value_col):
    """Return (labels, arrays) for the three fixed bins."""
    labels, arrays = [], []
    for label, grades in BINS:
        v = sub.loc[sub["grade"].isin(grades), value_col].dropna().values
        v = v[np.isfinite(v)]
        labels.append(label)
        arrays.append(v)
    return labels, arrays


stat_rows = []
raw_rows = []
panel_cache = {}

for var in VARIANTS:
    for tox, tox_title in TOXICITIES:
        sub = ici[(ici["toxicity"] == tox) & var["mask"](ici)].copy()
        labels, arrays = collect(sub, var["value"])

        adj_p = [None]
        for i in range(1, len(arrays)):
            a, b = arrays[i - 1], arrays[i]
            if len(a) >= 3 and len(b) >= 3:
                adj_p.append(mannwhitneyu(a, b, alternative="two-sided").pvalue)
            else:
                adj_p.append(None)

        jt_stat, jt_z, jt_p = jonckheere_terpstra(arrays)
        tau_b, tau_p, pct_tied = trend_tau(arrays)
        panel_cache[(var["key"], tox)] = (labels, arrays, adj_p, jt_z, jt_p,
                                          tau_b, tau_p, pct_tied)

        # Denominator for the "any post steroid" proportion is always the full
        # ICI+ group for that bin, independent of the variant's own mask.
        for i, (label, grades) in enumerate(BINS):
            all_ici = ici[(ici["toxicity"] == tox) & ici["grade"].isin(grades)]
            d = arrays[i]
            stat_rows.append({
                "variant": var["key"],
                "toxicity": tox,
                "grade_bin": label,
                "n_plotted": len(d),
                "n_ici_total": len(all_ici),
                "pct_plotted": 100.0 * len(d) / len(all_ici) if len(all_ici) else np.nan,
                "median": float(np.median(d)) if len(d) else np.nan,
                "q1": float(np.percentile(d, 25)) if len(d) else np.nan,
                "q3": float(np.percentile(d, 75)) if len(d) else np.nan,
                "post_mgdays_median": float(all_ici["post_steroid"].median()),
                "pre_mgdays_median": float(all_ici["pre_steroid"].median()),
                "pct_any_post_steroid": float(100.0 * (all_ici["post_steroid"] > 0).mean()),
                "pct_zero_both_windows": float(100.0 * (~all_ici["any_steroid"]).mean()),
                "p_vs_prev_bin": adj_p[i],
                # Primary trend result: tie-corrected.
                "trend_tau_b": tau_b,
                "trend_tau_p": tau_p,
                "pct_values_tied_at_zero": pct_tied,
                # Retained for continuity with the notebook; unreliable when
                # pct_values_tied_at_zero is large. See jonckheere_terpstra().
                "trend_JT_stat_uncorrected": jt_stat,
                "trend_z_uncorrected": jt_z,
                "trend_p_uncorrected": jt_p,
            })

for tox, _ in TOXICITIES:
    for g in sorted(ici["grade"].unique()):
        sub = ici[(ici["toxicity"] == tox) & (ici["grade"] == g)]
        if len(sub) == 0:
            continue
        raw_rows.append({
            "toxicity": tox,
            "grade": int(g),
            "n_ici": len(sub),
            "n_any_steroid": int(sub["any_steroid"].sum()),
            "n_post_gt0": int((sub["post_steroid"] > 0).sum()),
            "pct_any_post_steroid": float(100.0 * (sub["post_steroid"] > 0).mean()),
            "post_mgdays_median": float(sub["post_steroid"].median()),
            "pre_mgdays_median": float(sub["pre_steroid"].median()),
            "log2_ratio_median_zeros_excluded": float(sub["log2_ratio"].median()),
            "log2_ratio_median_zeros_included": float(sub["ratio_all"].median()),
        })

stats_df = pd.DataFrame(stat_rows)
stats_df.to_csv(OUT_DIR / "Steroid_Sensitivity_3variants_stats.csv", index=False)
pd.DataFrame(raw_rows).to_csv(OUT_DIR / "Steroid_Sensitivity_by_raw_grade.csv", index=False)
print(f"\nSaved stats CSVs to {OUT_DIR}")


# ---------------------------------------------------------------------------
# FIGURE
# ---------------------------------------------------------------------------
colors = bin_colors()
fig, axes = plt.subplots(3, 3, figsize=(9.2, 8.4), constrained_layout=True)

for r, var in enumerate(VARIANTS):
    for c, (tox, tox_title) in enumerate(TOXICITIES):
        ax = axes[r, c]
        (labels, arrays, adj_p, jt_z, jt_p,
         tau_b, tau_p, pct_tied) = panel_cache[(var["key"], tox)]
        positions = list(range(len(arrays)))

        nonempty = [(p, d, col) for p, d, col in zip(positions, arrays, colors) if len(d) > 1]
        if nonempty:
            parts = ax.violinplot([d for _, d, _ in nonempty],
                                  positions=[p for p, _, _ in nonempty],
                                  showmeans=False, showmedians=False,
                                  showextrema=False, widths=0.6)
            for pc, (_, _, col) in zip(parts["bodies"], nonempty):
                pc.set_facecolor(col)
                pc.set_alpha(0.35)
                pc.set_edgecolor(col)
                pc.set_linewidth(1.2)

        rng = np.random.default_rng(42)
        for pos, d, col in zip(positions, arrays, colors):
            if len(d) == 0:
                continue
            # Subsample the huge control group so the cloud stays legible.
            show = d if len(d) <= 2500 else rng.choice(d, 2500, replace=False)
            jitter = rng.normal(pos, 0.045, len(show))
            ax.scatter(jitter, show, alpha=0.35, s=5, color=col,
                       edgecolors="white", linewidths=0.15, zorder=3)
            ax.hlines(np.median(d), pos - 0.24, pos + 0.24,
                      colors="black", linewidth=1.7, zorder=4)

        allv = np.concatenate([d for d in arrays if len(d)])
        lo, hi = float(np.percentile(allv, 0.5)), float(np.percentile(allv, 99.5))
        if hi <= lo:
            hi = lo + 1.0
        span = hi - lo
        base = hi + span * 0.10
        tick = span * 0.025

        for i in range(1, len(arrays)):
            if adj_p[i] is None:
                continue
            x1, x2 = positions[i - 1] + 0.07, positions[i] - 0.07
            ax.plot([x1, x1, x2, x2], [base - tick, base, base, base - tick],
                    color="black", linewidth=0.9, zorder=5)
            ax.text((x1 + x2) / 2, base + span * 0.015, ptxt(adj_p[i]),
                    ha="center", va="bottom", fontsize=5.5, zorder=5)

        ax.set_xticks(positions)
        ax.set_xticklabels([f"{l}\nn={len(d):,}" for l, d in zip(labels, arrays)],
                           fontsize=6)
        ax.tick_params(axis="y", labelsize=6)
        ax.set_ylim(lo - span * 0.20, base + span * 0.20)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.axhline(0, color="black", linewidth=0.5, linestyle=":", alpha=0.5, zorder=1)

        trend = f"Trend tau$_b$={tau_b:.3f}, {ptxt(tau_p)}"
        ax.text(0.02, 1.005, trend, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=5.5, style="italic")
        if pct_tied > 20:
            ax.text(0.98, 1.005, f"{pct_tied:.0f}% tied at 0", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=5, style="italic", color="#b03a2e")

        if r == 0:
            ax.set_title(tox_title, fontsize=9, fontweight="bold", pad=14)
        if c == 0:
            ax.set_ylabel(var["ylabel"], fontsize=7)
            ax.text(-0.42, 0.5, var["row_label"], transform=ax.transAxes,
                    ha="center", va="center", fontsize=7.5, fontweight="bold",
                    rotation=90, linespacing=1.5)

fig.suptitle("Steroid exposure by irAE grade — three metric definitions, ICI+ patients",
             fontsize=10.5, fontweight="bold")

pdf_path = OUT_DIR / "Steroid_Sensitivity_3variants.pdf"
with PdfPages(str(pdf_path)) as pdf:
    pdf.savefig(fig, dpi=450)
fig.savefig(OUT_DIR / "Steroid_Sensitivity_3variants.png", dpi=300, bbox_inches="tight")
print(f"Saved figure: {pdf_path}")


# ---------------------------------------------------------------------------
# COMPANION FIGURE: the zeros used as a proportion rather than an imputed ratio
# ---------------------------------------------------------------------------
fig2, axes2 = plt.subplots(1, 3, figsize=(7.6, 2.4), constrained_layout=True)

for c, (tox, tox_title) in enumerate(TOXICITIES):
    ax = axes2[c]
    pcts, ns, ks = [], [], []
    for label, grades in BINS:
        grp = ici[(ici["toxicity"] == tox) & ici["grade"].isin(grades)]
        n = len(grp)
        k = int((grp["post_steroid"] > 0).sum())
        pcts.append(100.0 * k / n if n else np.nan)
        ns.append(n)
        ks.append(k)

    positions = list(range(len(BINS)))
    ax.bar(positions, pcts, width=0.62, color=colors, alpha=0.85,
           edgecolor="black", linewidth=0.6)

    # Wilson 95% interval — normal approximation is unreliable at these n.
    for pos, p, n in zip(positions, pcts, ns):
        if not n:
            continue
        ph = p / 100.0
        z = 1.96
        denom = 1 + z**2 / n
        centre = (ph + z**2 / (2 * n)) / denom
        half = z * np.sqrt(ph * (1 - ph) / n + z**2 / (4 * n**2)) / denom
        ax.plot([pos, pos], [100 * (centre - half), 100 * (centre + half)],
                color="black", linewidth=1.0)

    for pos, p, k, n in zip(positions, pcts, ks, ns):
        ax.text(pos, p + 3.5, f"{p:.0f}%", ha="center", va="bottom", fontsize=6.5,
                fontweight="bold")
        ax.text(pos, 3, f"{k:,}/{n:,}", ha="center", va="bottom", fontsize=5,
                color="white")

    ax.set_xticks(positions)
    ax.set_xticklabels([l for l, _ in BINS], fontsize=6.5)
    ax.set_ylim(0, 100)
    ax.tick_params(axis="y", labelsize=6)
    ax.set_title(tox_title, fontsize=8.5, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if c == 0:
        ax.set_ylabel("Patients with any steroid\nin post window (%)", fontsize=7)

fig2.suptitle("Proportion of ICI+ patients receiving any post-window steroid, by irAE grade",
              fontsize=9.5, fontweight="bold")

pdf2 = OUT_DIR / "Steroid_Sensitivity_proportion_treated.pdf"
with PdfPages(str(pdf2)) as pdf:
    pdf.savefig(fig2, dpi=450)
fig2.savefig(OUT_DIR / "Steroid_Sensitivity_proportion_treated.png", dpi=300,
             bbox_inches="tight")
print(f"Saved figure: {pdf2}")


# ---------------------------------------------------------------------------
# CONSOLE SUMMARY
# ---------------------------------------------------------------------------
print("\n" + "=" * 96)
print("TREND SUMMARY across No AE -> Grade 1-2 -> Grade 3-4")
print("tau_b is the primary result (tie-corrected). JT z is the notebook's")
print("uncorrected statistic and is invalid where %tied is large.")
print("=" * 96)
print(f"{'variant':<28}{'toxicity':<24}{'n':>8}{'%tied':>7}"
      f"{'tau_b':>8}{'tau p':>11}{'JT z':>8}")
for var in VARIANTS:
    for tox, _ in TOXICITIES:
        (labels, arrays, adj_p, jt_z, jt_p,
         tau_b, tau_p, pct_tied) = panel_cache[(var["key"], tox)]
        n = sum(len(d) for d in arrays)
        flag = "  <-- JT unreliable" if pct_tied > 20 else ""
        print(f"{var['key']:<28}{tox:<24}{n:>8,}{pct_tied:>6.0f}%"
              f"{tau_b:>8.3f}{tau_p:>11.2e}{jt_z:>8.2f}{flag}")

print("\n" + "=" * 78)
print("MEDIANS BY BIN")
print("=" * 78)
for var in VARIANTS:
    print(f"\n{var['key']}")
    for tox, _ in TOXICITIES:
        labels, arrays, adj_p = panel_cache[(var["key"], tox)][:3]
        cells = "  ".join(
            f"{l}: {np.median(d):6.2f} (n={len(d):,})" if len(d) else f"{l}: n/a"
            for l, d in zip(labels, arrays)
        )
        print(f"  {tox:<24}{cells}")

print("\n" + "=" * 78)
print("PROPORTION WITH ANY POST-WINDOW STEROID (the zeros, used as a proportion)")
print("=" * 78)
prop = (stats_df[stats_df["variant"] == "V1_ratio_zeros_excluded"]
        .pivot(index="toxicity", columns="grade_bin", values="pct_any_post_steroid"))
print(prop[["No AE", "Grade 1-2", "Grade 3-4"]].round(1).to_string())
