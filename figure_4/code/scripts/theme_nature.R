# theme_nature.R
# Reusable ggplot2 theme + helpers implementing the Nature-style publication
# figure rules documented in SKILL.md. Companion script for the
# `nature_publication_figures_r` skill.
#
# Usage:
#   source("theme_nature.R")
#   ggplot(df, aes(x, y, colour = group)) +
#     geom_point() +
#     scale_color_nature() +
#     theme_nature() +
#     nature_legend(n_entries = 3)

suppressPackageStartupMessages(library(ggplot2))

# Nature rule: prefer Arial (or Helvetica). We register Arial so its family name
# resolves consistently, and default `theme_nature(base_family = "Arial")` to it.
# On macOS Arial ships as a system font; on Linux install "ttf-mscorefonts" or
# substitute Helvetica/Nimbus Sans and pass that as base_family.
if (requireNamespace("systemfonts", quietly = TRUE)) {
  .arial <- tryCatch({
    if ("match_fonts" %in% getNamespaceExports("systemfonts")) {
      systemfonts::match_fonts("Arial")$path          # systemfonts >= 1.1.0
    } else {
      systemfonts::match_font("Arial")$path           # older systemfonts
    }
  }, error = function(e) NULL)
  if (!is.null(.arial) && nzchar(.arial)) {
    try(systemfonts::register_font(name = "Arial", plain = .arial), silent = TRUE)
  }
}

# Okabe-Ito colorblind-safe palette (accessible under CVD, strong contrast).
nature_palette <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

scale_color_nature <- function(...) ggplot2::scale_color_manual(values = nature_palette, ...)
scale_fill_nature  <- function(...) ggplot2::scale_fill_manual(values = nature_palette, ...)

# ggplot2 line widths are in mm; `.pt` converts pt -> ggplot's internal unit.
pt_to_lwd <- 1 / .pt

#' Nature-style ggplot2 theme
#'
#' @param base_size base text size in pt for body text (rule: 5-6 pt)
#' @param base_family font family; must be installed and registered
#'   (e.g. via showtext/systemfonts) — rule prefers Arial or Helvetica
#' @param panel_letter_size size in pt for bold panel/strip labels (rule: 8 pt)
theme_nature <- function(base_size = 6, base_family = "Arial", panel_letter_size = 8) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = pt_to_lwd),
      axis.ticks       = element_line(colour = "black", linewidth = pt_to_lwd),
      # Negative length pushes ticks outward, since ggplot2 has no native option.
      axis.ticks.length = unit(-2.75, "pt"),
      axis.text.x  = element_text(margin = margin(t = 6), colour = "black", size = base_size),
      axis.text.y  = element_text(margin = margin(r = 6), colour = "black", size = base_size),
      axis.title   = element_text(size = base_size, colour = "black"),
      text         = element_text(colour = "black"),
      plot.title   = element_blank(),
      strip.background = element_blank(),
      strip.text   = element_text(size = panel_letter_size, face = "bold", hjust = 0),
      legend.background = element_blank(),
      legend.key   = element_blank(),
      legend.text  = element_text(size = base_size),
      legend.title = element_text(size = base_size)
    )
}

#' Legend placement per Nature rules: above (1 row) if <5 entries, else right (1 column).
nature_legend <- function(n_entries) {
  if (n_entries < 5) {
    theme(legend.position = "top", legend.direction = "horizontal")
  } else {
    theme(legend.position = "right", legend.direction = "vertical")
  }
}

#' Figure size in Nature's 1.5-inch layout units, warns past the 5 x 7 in caps.
#'
#' @param width_units number of 1.5in width units
#' @param height_units number of 1.5in height units
#' @return list(width, height) in inches, for use with ggsave()
nature_size <- function(width_units = 2, height_units = 2) {
  w <- width_units * 1.5
  h <- height_units * 1.5
  if (w > 5) warning(sprintf("Total width %.2gin exceeds the 5in Nature limit.", w))
  if (h > 7) warning(sprintf("Total height %.2gin exceeds the 7in Nature limit.", h))
  list(width = w, height = h)
}

#' FDR-correct a table of test results and attach Nature-style significance labels.
#'
#' Choose the underlying test yourself (paired/unpaired, parametric/non-parametric,
#' 2 groups vs >2 groups) — this only pools the p-values for one shared correction
#' and formats compact labels for annotation.
#'
#' @param test_results data.frame with a raw p-value column
#' @param p_col name of the raw p-value column (default "p")
#' @return test_results with `p.adj` and `label` columns added
fdr_annotate <- function(test_results, p_col = "p") {
  test_results$p.adj <- p.adjust(test_results[[p_col]], method = "BH")
  test_results$label <- cut(
    test_results$p.adj,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
  )
  test_results
}

#' Is a working cairo device available? (needs XQuartz on macOS)
#'
#' Opening a cairo_pdf can succeed on a machine where cairo still cannot resolve
#' a real font — the failure ("failed to find or load PDF CID font") only
#' surfaces once text is measured. So actually measure text in `family` here,
#' otherwise save_nature_figure() picks cairo and then dies at render time.
.cairo_ok <- function(family = "Arial") {
  isTRUE(tryCatch({
    f <- tempfile(fileext = ".pdf")
    grDevices::cairo_pdf(f)
    on.exit({
      try(grDevices::dev.off(), silent = TRUE)
      unlink(f)
    }, add = TRUE)
    graphics::par(family = family)
    w <- graphics::strwidth("Ag", units = "inches")
    is.finite(w) && w > 0
  }, error = function(e) FALSE, warning = function(w) FALSE))
}

#' Pick a vector device that produces editable RGB text with system fonts.
#'
#' Order of preference:
#'   1. cairo_pdf / cairo_ps        (cross-platform, needs cairo/XQuartz)
#'   2. macOS quartz(type="pdf")    (native Core Graphics, no XQuartz needed)
#'   3. base pdf()                  (last resort; Arial embedding may fall back)
#' All three keep text as editable text and render RGB.
.nature_pdf_device <- function(eps = FALSE) {
  if (.cairo_ok()) {
    if (eps) grDevices::cairo_ps else grDevices::cairo_pdf
  } else if (Sys.info()[["sysname"]] == "Darwin") {
    function(filename, ...) grDevices::quartz(file = filename, type = "pdf", ...)
  } else {
    function(filename, ...) grDevices::pdf(file = filename, useDingbats = FALSE, ...)
  }
}

#' Save a figure per Nature export rules: vector PDF/EPS, RGB, editable text,
#' >=450 dpi for any raster content.
#'
#' @param plot a ggplot/patchwork object
#' @param filename output path; should end in .pdf or .eps
#' @param width_units,height_units passed to nature_size()
#' @param dpi raster resolution for any embedded raster layers (rule: >=450)
save_nature_figure <- function(plot, filename, width_units = 2, height_units = 2, dpi = 450) {
  sz <- nature_size(width_units, height_units)
  is_eps <- grepl("\\.eps$", filename, ignore.case = TRUE)
  if (!grepl("\\.(pdf|eps)$", filename, ignore.case = TRUE)) {
    warning("Nature figures should export as vector PDF/EPS; got: ", filename)
  }
  ggsave(
    filename, plot,
    width = sz$width, height = sz$height, units = "in",
    device = .nature_pdf_device(eps = is_eps),
    dpi = dpi
  )
}
