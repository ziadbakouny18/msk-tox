# Code for "Population-Scale Precision Safety in Oncology Reveals Clinical and Genetic Determinants of Systemic Therapy Toxicity"

## System requirements

### Operating systems
- Linux x86_64 (analyses run on a SLURM HPC cluster)
- macOS 14+ (figure generation)

### Software
- **Python 3.14.7:** pandas, numpy, scipy, matplotlib, seaborn, scikit-learn, statsmodels, lifelines, tqdm, Pillow
  - PyMOL (structural rendering, Figure 5d only)
- **R 4.6.1:** tidyverse (dplyr, tidyr, readr, tibble, stringr, lubridate), ggplot2, patchwork, ggrepel, scales, systemfonts, svglite, ragg, data.table, jsonlite, httr2, survival, cmprsk, coxme
  - SPACox is not on CRAN; install it from [github.com/WenjianBI/SPACox](https://github.com/WenjianBI/SPACox)
- **Command-line tools:** PLINK 1.9, PLINK 2.0, bcftools, curl

### Non-standard hardware
- SLURM-managed HPC cluster (genome-wide association analyses)
- GPU (HLA allele imputation)

## Repository layout

Each figure has its own top-level directory with a README describing its panels
and expected inputs. No source data is included except where noted.

| Directory | Contents |
|---|---|
| `figure_1/code/` | RAG-LLM validation, laboratory and corticosteroid analyses |
| `figure_2/code/` | Toxicity incidence, treatment/demographic and outcome analyses |
| `figure_3/plotting/` | Random survival forest models, rendered from `figure_3/source_data/` |
| `figure_4/code/scripts/` | Genome-wide association pipeline, replication and annotation |
| `figure_5/` | HLA-DRB1 positional analyses and structural visualization |

`figure_3` is the one figure that ships source data: aggregate CSVs in
`source_data/` plus a `MANIFEST.csv`, with renderers in `plotting/`. The Figure 4
scripts are nested one level deeper, in `code/scripts/`, because several of them
resolve their own helpers through a literal `scripts/` path component.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0). See the [LICENSE](LICENSE) file for the full text.

Methods described herein are the subject of a pending patent application. For licensing inquiries, contact Memorial Sloan Kettering Cancer Center at [tech transfer email].
