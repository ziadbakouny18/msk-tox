# PREP: turn the raw FinnGen per-variant PheWAS JSONs into one tidy table.
# Output: results/figS4/finngen_phewas.tsv (consumed by phewas_plot.R)

library(jsonlite)
library(tidyverse)

root <- normalizePath(".", mustWork = TRUE)
fig <- file.path(root, "results", "figS4")

files <- list.files(file.path(fig, "finngen"), pattern = "_finngen\\.json$", full.names = TRUE)

ph <- map_dfr(files, function(f) {
  rs <- str_extract(basename(f), "^[^_]+")
  d  <- fromJSON(f, simplifyDataFrame = TRUE)
  r  <- d$results
  if (is.null(r) || length(r) == 0) return(tibble())
  as_tibble(r) %>%
    transmute(rsid = rs, phenocode, phenostring, category,
              pval, beta, sebeta, mlogp, n_case, n_control)
}) %>%
  filter(!is.na(pval))

write_tsv(ph, file.path(fig, "finngen_phewas.tsv"))
message("wrote finngen_phewas.tsv: ", nrow(ph), " endpoint associations across ",
        n_distinct(ph$rsid), " variants")
