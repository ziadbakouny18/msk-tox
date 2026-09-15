# figure_3 — MSK-Tox reproducible figures

This folder reproduces Figure 3 and Supplementary Figure 3 of the MSK-Tox
paper from source data.

## Layout

| path | holds |
|---|---|
| `source_data/` | source CSVs with aggregate statistics plus `MANIFEST.csv` |
| `plotting/` | 4 self-contained renderers plus vendored `aesthetics.py` |
| `figures/` | Rendered PNG and PDF outputs |

## Setup

Requires Python 3.12

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run every command from this directory (`msk-tox-paper/figure_3`). Each panel
command writes both PNG and PDF into `figures/`.

## Replication recipes

**fig3b** — ICI vs non-ICI summary over the five `_grade0_base` models
(adrenal insufficiency, hypothyroidism, liver toxicity, pneumonitis, colitis).
Reads `ici_nonici_cindex.csv`, `ici_nonici_incidence_curves.csv`,
`ici_nonici_logrank.csv`, `ici_nonici_panel_meta.csv`.

```bash
python plotting/plot_ici_nonici_comparison.py --panel_set A
```

**fig3c** — same summary on the pan-AE models. Reads the same four CSVs.

```bash
python plotting/plot_ici_nonici_comparison.py --panel_set B   # all-AE grade 0
python plotting/plot_ici_nonici_comparison.py --panel_set C   # all-AE grade 3
```

**supp3a** — stratified c-index barplot over the five `_grade0_base` models,
stratified by line of therapy. Reads `cindex_stratified_lot.csv`.

```bash
python plotting/plot_cindex_barplots_stratified.py --stratify_by lot
```

**supp3b** — same, stratified by treatment modality. Reads
`cindex_stratified_treatment.csv`.

```bash
python plotting/plot_cindex_barplots_stratified.py --stratify_by treatment
```

**supp3c** — SHAP beeswarm for liver toxicity and hypothyroidism in Renal
Cell Carcinoma. NOTE: Data is not provided here as patient-level data is required to reproduce the plot.

```bash
python plotting/plot_shap_beeswarm.py --ae liver_toxicity
python plotting/plot_shap_beeswarm.py --ae hypothyroidism
```

**supp3e** — c-index swarm over all `_grade0` models across base, +labs,
+genomics, +labs+genomics. Reads `cindex_swarm.csv`.

```bash
python plotting/plot_cindex_swarm.py
```
