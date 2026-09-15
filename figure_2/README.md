# figure_2 — MSK-Tox Figure 2 and Supplementary Figure 2

Scripts that generate Figure 2 (2A–2F) and Supplementary Figure 2 (S2A–S2G).

No source data is included here. Inputs are expected at
`figures/figures_data/figure 2/data` (same layout as the working OneDrive tree).
Notebooks assume they are run from a `figure 2/scripts/` equivalent of that tree.
Panels are written to `figure 2/results/main` and `figure 2/results/supp/<panel>`.

## Layout

| path | holds |
|---|---|
| `code/` | Figure 2 and S2 notebooks, plus the R frailty-model script |

## Panels

**Figure 2:** `Cumulative_Incidence_2A.ipynb`, `One_Year_Cumulative_Incidence_2B.ipynb`, `Cumulative_Incidence_By_ICI_&_Treatment_2C_D.ipynb` (produces both 2C and 2D), `Heatmap_CI_By_Cancer_Type_2E.ipynb`, `Forest_Cox_PFS_2F.ipynb`.

**Supplementary Figure 2:** `Cancer_Type_Treatment_S2A.ipynb`, `Grade_Histogram_S2B.ipynb`, `Sex_Prevalance_S2C.ipynb`, `Age_Group_Prevalance_LLM_S2D.ipynb`, `LOT_Prevalence_S2E.ipynb`, `Forest_Cox_OS_S2F.ipynb`, `Grade_Distribution_G3__S2G.ipynb` (plus `lot_frailty_cox_models.R`).

`lot_frailty_cox_models.R` fits the shared-frailty and clustered-robust Cox models for
S2E. Run `LOT_Prevalence_S2E.ipynb` first — it exports the per-line-of-therapy table the
R script reads. Run the R script from the `scripts/` folder so its relative paths resolve.
