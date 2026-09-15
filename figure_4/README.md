# figure_4 — MSK-Tox Figure 4 and Supplementary Figure 4

Scripts that generate Figure 4 and Supplementary Figure 4 (S4): the genome-wide
association pipeline, DFCI replication, LD-proxy lookups, and the eQTL, ancestry,
and PheWAS panels.

No source data is included here. Inputs are expected in the analysis trees the
scripts were run against (MSK OneDrive, Box, and the SLURM project directory), so
paths must be repointed before rerunning.

Run the scripts from `code/`. Several of them resolve helpers through a literal
`scripts/` path component — for example `source(file.path("scripts", "theme_nature.R"))`
— so the inner `scripts/` folder is kept as-is and the scripts are unmodified.

## Layout

| path | holds |
|---|---|
| `code/` | Figure 4 and S4 scripts, as the `scripts/` package tree they expect |

Within `code/scripts/`:

| path | holds |
|---|---|
| `00_gwas_pipeline/` | SPACox AE GWAS prep, worker, and SLURM array submission |
| `02_fig4_forest/` | Figure 4 competing-risk cumulative incidence and by-cohort forest plots |
| `03_dfci_replication/` | Streaming QC scan over the DFCI (PROFILE) full GWAS results |
| `05_ld_proxies/` | 1000G LD-proxy finders and the lead-SNP list |
| `08_figS4_eqtl_ancestry/` | S4 GTEx eQTL and 1000G ancestry panels, plus shared `R/panel_helpers.R` |
| `09_figS4_phewas/` | S4 FinnGen PheWAS prep and plot |

## Panels

**Figure 4:** `02_fig4_forest/snp_ci_forest.R`.

**Supplementary Figure 4:** `08_figS4_eqtl_ancestry/panel_foxe1_gtex.R`,
`panel_foxe1_ancestry.R`, `panel_foxp2_left_ventricle_gtex.R`,
`09_figS4_phewas/phewas_prep.R` and `phewas_plot.R`.

**Supplementary tables:** `01_dfci_gwas_replication_table.R` (DFCI replication) and
`02_prior_loci_ld_proxy_table.R` (prior pooled-irAE loci, Groha et al. 2022).

`theme_nature.R` is the shared ggplot2 theme sourced by the panel scripts.

## Notes

`00_gwas_pipeline` runs on a SLURM cluster and requires SPACox, which is not on
CRAN — install it from [github.com/WenjianBI/SPACox](https://github.com/WenjianBI/SPACox).

`05_ld_proxies/run_ld_proxies.sh` streams 1000G Phase 3 (hg19) over HTTPS and
computes proxies with PLINK 1.9, keeping all intermediates inside its own directory.

Lead variants analysed: rs112463084 (chr6), rs112386546 (chr8, *ESRP1*),
rs10983761 (chr9, *FOXE1*), rs12705907 (chr7, *FOXP2*).
