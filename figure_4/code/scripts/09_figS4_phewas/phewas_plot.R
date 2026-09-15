#!/usr/bin/env Rscript

# Publication-quality FinnGen PheWAS panel for rs10983761 (FOXE1).
# Run from the manuscript_reproduction_2026-08-23 package root.
# The four strongest thyroid endpoints whose full names fit are labelled.

# ---- Packages ----

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ragg)
})

source(file.path("scripts", "theme_nature.R"))

# ---- Settings ----

root <- normalizePath(".", mustWork = TRUE)
fig_dir <- file.path(root, "results", "figS4")
panel_dir <- file.path(fig_dir, "panels")

RS <- "rs10983761"
GENE <- "FOXE1"
HILITE_RE <- "thyroid|hypothyroid|goit|thyreo|graves"
HILITE_LAB <- "Thyroid disease"
N_LABEL <- 4
MAXCHAR <- 46
FONT <- "Arial"

dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data ----

# Map FinnGen registry groupings to compact phenotype chapters.
catname <- function(s) {
  code <- ifelse(grepl("\\(", s), sub(".*\\(([^)]*)\\)\\s*$", "\\1", s), "")

  case_when(
    grepl("AB1", code) ~ "Infections",
    grepl("CD2|ICD-O-3", code) | grepl("Neoplasm", s) ~ "Neoplasms",
    grepl("^D3", code) ~ "Blood/immune",
    grepl("^E4", code) ~ "Endocrine/metabolic",
    grepl("Diabetes", s) ~ "Diabetes",
    grepl("^F5", code) ~ "Psychiatric",
    grepl("^G6", code) ~ "Neurologic",
    grepl("^H7", code) ~ "Eye",
    grepl("^H8", code) ~ "Ear",
    grepl("^I9", code) ~ "Cardiovascular",
    grepl("^J10", code) ~ "Respiratory",
    grepl("^K11", code) ~ "Digestive",
    grepl("Gastrointestinal", s) ~ "Gastrointestinal",
    grepl("^L12", code) ~ "Skin",
    grepl("^M13", code) ~ "Musculoskeletal",
    grepl("Rheuma", s) ~ "Rheuma",
    grepl("^N14", code) ~ "Genitourinary",
    grepl("^O15", code) ~ "Pregnancy",
    grepl("^P16", code) ~ "Perinatal",
    grepl("^Q17", code) ~ "Congenital",
    TRUE ~ "Other"
  )
}

cat_order <- c(
  "Infections", "Neoplasms", "Blood/immune", "Endocrine/metabolic",
  "Diabetes", "Psychiatric", "Neurologic", "Eye", "Ear",
  "Cardiovascular", "Respiratory", "Digestive", "Gastrointestinal",
  "Skin", "Musculoskeletal", "Rheuma", "Genitourinary", "Pregnancy",
  "Perinatal", "Congenital", "Other"
)

ph <- read_tsv(
  file.path(fig_dir, "finngen_phewas.tsv"),
  show_col_types = FALSE
)

stopifnot(all(c(
  "rsid", "phenocode", "phenostring", "category", "pval", "beta",
  "sebeta", "n_case", "n_control"
) %in% names(ph)))

ph <- ph |>
  filter(rsid == RS, !is.na(pval), pval > 0) |>
  mutate(cat = factor(catname(category), levels = cat_order)) |>
  arrange(cat, desc(pval)) |>
  mutate(
    x = row_number(),
    neglog10p = -log10(pval),
    group = if_else(
      grepl(HILITE_RE, paste(phenostring, category), ignore.case = TRUE),
      HILITE_LAB,
      "Other phenotype"
    )
  )

stopifnot(nrow(ph) > 0, all(is.finite(ph$neglog10p)))
glimpse(ph)

# Alternate neutral tones by adjacent chapter; reserve orange and triangles for
# thyroid endpoints so the highlight remains legible without relying on colour.
present <- ph |>
  distinct(cat) |>
  arrange(cat) |>
  mutate(chapter_tone = if_else(row_number() %% 2 == 0, "greyA", "greyB"))

ph <- ph |>
  left_join(present, by = "cat") |>
  mutate(
    colour_group = if_else(group == HILITE_LAB, HILITE_LAB, chapter_tone),
    shape_group = if_else(group == HILITE_LAB, HILITE_LAB, "Other phenotype")
  )

# ---- Labels ----

# Remove registry boilerplate while preserving complete endpoint names.
clean_name <- function(s) {
  s <- gsub(
    "\\s*,?\\s*\\(?(excluding all cancers|controls excluding all cancers)\\)?",
    "",
    s,
    ignore.case = TRUE
  )
  s <- trimws(gsub("\\s+", " ", s))
  sub("[,;:]+$", "", s)
}

lab <- ph |>
  filter(group == HILITE_LAB) |>
  arrange(pval) |>
  mutate(nm = clean_name(phenostring)) |>
  filter(nchar(nm) <= MAXCHAR) |>
  slice_head(n = N_LABEL) |>
  mutate(ltxt = paste0(nm, " (P = ", formatC(pval, format = "e", digits = 0), ")"))

stopifnot(nrow(lab) == N_LABEL)

ph <- ph |>
  mutate(labeled = phenocode %in% lab$phenocode)

mids <- ph |>
  group_by(cat) |>
  summarise(mid = mean(x), .groups = "drop")

# ---- Plot ----

y_max <- max(ph$neglog10p) * 1.14

p <- ggplot(ph, aes(x = x, y = neglog10p)) +
  geom_point(
    aes(colour = colour_group, shape = shape_group),
    size = 1.05,
    alpha = 0.9,
    stroke = 0.15
  ) +
  geom_text_repel(
    data = lab,
    aes(label = ltxt),
    family = FONT,
    colour = "black",
    size = 1.75,
    seed = 17,
    box.padding = 0.55,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.colour = "grey45",
    force = 5,
    force_pull = 0.35,
    max.iter = 50000,
    max.time = 3,
    max.overlaps = Inf,
    direction = "both"
  ) +
  scale_colour_manual(
    values = c(
      "Thyroid disease" = "#D55E00",
      "greyA" = "#BDBDBD",
      "greyB" = "#4D4D4D"
    ),
    breaks = HILITE_LAB,
    labels = HILITE_LAB,
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Thyroid disease" = 17, "Other phenotype" = 16),
    breaks = HILITE_LAB,
    labels = HILITE_LAB,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = mids$mid,
    labels = mids$cat,
    expand = expansion(add = 3),
    guide = guide_axis(angle = 90)
  ) +
  scale_y_continuous(
    limits = c(0, y_max),
    expand = expansion(mult = c(0, 0.01)),
    breaks = scales::breaks_pretty(n = 4)
  ) +
  labs(
    x = "Phenotype chapter",
    y = expression(-log[10] ~ italic(P))
  ) +
  guides(
    colour = guide_legend(override.aes = list(shape = 17, size = 2, alpha = 1)),
    shape = "none"
  ) +
  theme_nature(base_size = 6, base_family = FONT) +
  nature_legend(n_entries = 1) +
  theme(
    axis.text.x = element_text(size = 5.2, hjust = 1, vjust = 0.5),
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(9, "pt"),
    plot.margin = margin(4, 5, 4, 5)
  )

# Show the plot in RStudio without opening an implicit PDF device under Rscript.
if (interactive()) {
  print(p)
}

# ---- Export ----

# Preserve the original manuscript panel dimensions requested for this figure.
# This intentionally exceeds the generic 5-inch single-panel Nature guideline.
size <- list(width = 8, height = 3)

pdf_file <- file.path(panel_dir, paste0("phewas_", RS, ".pdf"))
svg_file <- file.path(panel_dir, paste0("phewas_", RS, ".svg"))
png_file <- file.path(panel_dir, paste0("phewas_", RS, ".png"))

save_nature_figure(
  p,
  pdf_file,
  width_units = 16 / 3,
  height_units = 2,
  dpi = 600
)

ggsave(
  svg_file,
  p,
  width = size$width,
  height = size$height,
  units = "in",
  device = svglite::svglite
)

ggsave(
  png_file,
  p,
  width = size$width,
  height = size$height,
  units = "in",
  dpi = 600,
  device = ragg::agg_png,
  background = "white"
)

ph |>
  transmute(
    rsid,
    phenocode,
    phenostring,
    category,
    chapter = cat,
    pval,
    beta,
    sebeta,
    n_case,
    n_control,
    neglog10p,
    group,
    labeled
  ) |>
  arrange(pval) |>
  write_csv(file.path(fig_dir, paste0("phewas_", RS, "_plotdata.csv")))

message(
  "Wrote publication PheWAS panel: ",
  paste(basename(c(pdf_file, svg_file, png_file)), collapse = ", "),
  " | endpoints=", nrow(ph),
  " | labels=", sum(ph$labeled)
)
