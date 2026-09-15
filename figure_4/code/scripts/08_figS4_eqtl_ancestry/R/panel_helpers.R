suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(httr2)
  library(jsonlite)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# Run from the manuscript_reproduction_2026-08-23 package root. All cached
# public inputs and outputs remain inside this package.
root <- normalizePath(".", mustWork = TRUE)
out_root <- file.path(root, "results", "figS4")
cache_dir <- file.path(out_root, "gtex_cache")
table_dir <- out_root
plot_dir <- file.path(out_root, "panels")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

first_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[[1]]
}

lead_snps <- read_tsv(
  file.path(root, "scripts", "05_ld_proxies", "lead_snps.tsv"),
  show_col_types = FALSE
)

# I only keep the panels that are useful for the final figure.
panel_targets <- tribble(
  ~lead_snp,      ~trait,                   ~target_gene, ~gencode_id,            ~gtex_tissue,
  "rs10983761",   "Hypothyroidism",         "FOXE1",      "ENSG00000178919.8",  "Thyroid",
  "rs12705907",   "Pneumonitis",            "FOXP2",      "ENSG00000128573.24", "Lung",
  "rs12705907",   "Pneumonitis",            "FOXP2",      "ENSG00000128573.24", "Heart_Left_Ventricle",
  "rs112463084",  "Adrenal insufficiency",  "HLA-DRB1",   "ENSG00000196126.11", NA_character_
) %>%
  left_join(lead_snps, by = c("lead_snp" = "snp"))

target_row <- function(lead_snp, target_gene, gtex_tissue = NA_character_) {
  row <- panel_targets %>%
    filter(
      lead_snp == .env$lead_snp,
      target_gene == .env$target_gene,
      if (is.na(.env$gtex_tissue)) is.na(gtex_tissue) else gtex_tissue == .env$gtex_tissue
    )
  stopifnot(nrow(row) == 1)
  row[1, ]
}

cache_key <- function(...) {
  str_replace_all(paste(..., sep = "__"), "[^A-Za-z0-9_.-]+", "_")
}

fetch_gtex_json <- function(path, params, key) {
  cache_file <- file.path(cache_dir, paste0(key, ".json"))
  if (file.exists(cache_file) && file.size(cache_file) > 0) {
    cached <- fromJSON(cache_file, simplifyVector = FALSE)
    if (is.null(cached$detail) || !str_detect(cached$detail, "Could not resolve|Failed to perform HTTP request|HTTP 400")) {
      return(cached)
    }
  }

  req <- request(paste0("https://gtexportal.org", path)) %>%
    req_url_query(!!!params) %>%
    req_user_agent("ae-gwas-missing-panel-script/1.0") %>%
    req_error(is_error = function(resp) FALSE)

  resp <- tryCatch(req_perform(req), error = function(e) e)
  if (inherits(resp, "error")) {
    body <- conditionMessage(resp)
    writeLines(toJSON(list(detail = body), auto_unbox = TRUE), cache_file)
    return(list(detail = body))
  }

  body <- resp_body_string(resp)
  writeLines(body, cache_file)
  fromJSON(body, simplifyVector = FALSE)
}

extract_dyneqtl <- function(result, target) {
  if (!is.null(result$detail)) {
    return(tibble(
      lead_snp = target$lead_snp,
      trait = target$trait,
      target_gene = target$target_gene,
      tissue = target$gtex_tissue,
      gencode_id = target$gencode_id,
      available = FALSE,
      detail = result$detail
    ))
  }

  tibble(
    lead_snp = target$lead_snp,
    trait = target$trait,
    target_gene = result$geneSymbol %||% target$target_gene,
    tissue = result$tissueSiteDetailId %||% target$gtex_tissue,
    gencode_id = result$gencodeId %||% target$gencode_id,
    variant_id = result$variantId %||% NA_character_,
    available = TRUE,
    p_value = result$pValue %||% NA_real_,
    nes = result$nes %||% NA_real_,
    maf = result$maf %||% NA_real_,
    n_hom_ref = result$homoRefCount %||% NA_integer_,
    n_het = result$hetCount %||% NA_integer_,
    n_hom_alt = result$homoAltCount %||% NA_integer_,
    genotype = unlist(result$genotypes),
    expression = unlist(result$data)
  )
}

parse_gtex_variant <- function(variant_id) {
  parsed <- str_match(variant_id, "^chr[^_]+_[0-9]+_([^_]+)_([^_]+)_b[0-9]+$")
  tibble(ref = parsed[, 2], alt = parsed[, 3])
}

load_gtex_target <- function(target) {
  key <- cache_key("gtex_dyneqtl", "v8", target$lead_snp, target$target_gene, target$gtex_tissue)
  result <- fetch_gtex_json(
    "/api/v2/association/dyneqtl",
    list(
      tissueSiteDetailId = target$gtex_tissue,
      gencodeId = target$gencode_id,
      variantId = target$lead_snp,
      datasetId = "gtex_v8"
    ),
    key
  )

  extracted <- extract_dyneqtl(result, target)
  summary <- extracted %>%
    mutate(
      variant_id = if (!"variant_id" %in% names(.)) NA_character_ else variant_id,
      p_value = if (!"p_value" %in% names(.)) NA_real_ else p_value,
      nes = if (!"nes" %in% names(.)) NA_real_ else nes,
      maf = if (!"maf" %in% names(.)) NA_real_ else maf,
      n_hom_ref = if (!"n_hom_ref" %in% names(.)) NA_integer_ else n_hom_ref,
      n_het = if (!"n_het" %in% names(.)) NA_integer_ else n_het,
      n_hom_alt = if (!"n_hom_alt" %in% names(.)) NA_integer_ else n_hom_alt,
      detail = if (!"detail" %in% names(.)) NA_character_ else detail
    ) %>%
    group_by(lead_snp, trait, target_gene, tissue, gencode_id) %>%
    summarise(
      available = any(available),
      variant_id = first_or_na(variant_id),
      p_value = first_or_na(p_value),
      nes = first_or_na(nes),
      maf = first_or_na(maf),
      n_hom_ref = first_or_na(n_hom_ref),
      n_het = first_or_na(n_het),
      n_hom_alt = first_or_na(n_hom_alt),
      detail = first_or_na(detail),
      .groups = "drop"
    )

  points <- extracted %>%
    filter(available, !is.na(genotype), !is.na(expression)) %>%
    bind_cols(parse_gtex_variant(.$variant_id)) %>%
    mutate(
      genotype_dosage = as.integer(genotype),
      genotype_text = case_when(
        genotype_dosage == 0 ~ paste0(ref, ref),
        genotype_dosage == 1 ~ paste0(ref, alt),
        genotype_dosage == 2 ~ paste0(alt, alt),
        TRUE ~ NA_character_
      ),
      dosage_group = factor(genotype_dosage, levels = c(0, 1, 2))
    )

  list(summary = summary, points = points)
}

format_pairwise_p <- function(x) {
  case_when(
    is.na(x) ~ "p=NA",
    x < 0.001 ~ "p<0.001",
    TRUE ~ paste0("p=", signif(x, 2))
  )
}

format_nes <- function(x) {
  ifelse(is.na(x), "NES = NA", paste0("NES = ", formatC(x, format = "f", digits = 3)))
}

pairwise_gtex_tests <- function(points, genotype_counts, y_base, y_span) {
  if (nrow(genotype_counts) < 2) {
    return(tibble())
  }

  pairs <- t(combn(seq_len(nrow(genotype_counts)), 2)) %>%
    as_tibble(.name_repair = ~ c("x1", "x2")) %>%
    mutate(distance = x2 - x1) %>%
    arrange(distance, x1) %>%
    mutate(bracket_index = row_number())

  pairs %>%
    rowwise() %>%
    mutate(
      p_value = {
        g1 <- genotype_counts$genotype_x[x1]
        g2 <- genotype_counts$genotype_x[x2]
        y1 <- points$expression[points$genotype_x == g1]
        y2 <- points$expression[points$genotype_x == g2]
        if (length(y1) < 2 || length(y2) < 2) NA_real_ else wilcox.test(y1, y2)$p.value
      },
      y = y_base + (bracket_index - 1) * 0.115 * y_span,
      y_tip = y - 0.035 * y_span,
      label = format_pairwise_p(p_value)
    ) %>%
    ungroup()
}

make_gtex_plot <- function(target) {
  loaded <- load_gtex_target(target)
  summary <- loaded$summary
  points <- loaded$points
  title <- paste0(target$lead_snp, " / ", target$target_gene)
  # Tissue is named in the figure legend rather than on the panel itself, so the
  # subtitle carries only the effect size. The tissue is still recorded in the
  # output filename and in the accompanying *_gtex_summary.tsv.
  subtitle <- format_nes(first_or_na(summary$nes))

  if (nrow(points) == 0) {
    return(
      ggplot(tibble(x = 1, y = 1), aes(x, y)) +
        annotate("text", x = 1, y = 1, label = "Not available", size = 2.4) +
        labs(title = title, subtitle = subtitle) +
        theme_void(base_size = 7.5) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 8.2),
          plot.subtitle = element_text(hjust = 0.5, size = 7)
        )
    )
  }

  genotype_counts <- points %>%
    count(genotype_dosage, genotype_text, dosage_group, name = "n") %>%
    arrange(genotype_dosage) %>%
    mutate(
      genotype_x = factor(genotype_text, levels = genotype_text),
      x_label = paste0(genotype_text, "\nn=", comma(n))
    )

  points <- points %>%
    left_join(select(genotype_counts, genotype_dosage, genotype_x), by = "genotype_dosage")

  y_span <- diff(range(points$expression, na.rm = TRUE))
  if (y_span == 0) y_span <- 1
  bracket_df <- pairwise_gtex_tests(
    points,
    genotype_counts,
    y_base = max(points$expression, na.rm = TRUE) + 0.12 * y_span,
    y_span = y_span
  )
  y_upper <- if (nrow(bracket_df) > 0) max(bracket_df$y, na.rm = TRUE) + 0.12 * y_span else max(points$expression, na.rm = TRUE) + 0.12 * y_span

  p <- ggplot(points, aes(x = genotype_x, y = expression, fill = dosage_group, color = dosage_group)) +
    geom_violin(width = 0.82, trim = FALSE, alpha = 0.35, linewidth = 0.45) +
    geom_point(position = position_jitter(width = 0.11, height = 0, seed = 7), alpha = 0.42, size = 0.55) +
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.45, color = "black", linewidth = 0.34) +
    scale_x_discrete(labels = setNames(genotype_counts$x_label, genotype_counts$genotype_text)) +
    scale_fill_manual(values = c(`0` = "#b8c4c7", `1` = "#69b3e7", `2` = "#52d38a"), guide = "none") +
    scale_color_manual(values = c(`0` = "#b8c4c7", `1` = "#69b3e7", `2` = "#52d38a"), guide = "none") +
    labs(title = title, subtitle = subtitle, x = "Genotype", y = "Normalized expression") +
    coord_cartesian(ylim = c(min(points$expression, na.rm = TRUE) - 0.05 * y_span, y_upper), clip = "off") +
    theme_classic(base_size = 7.5) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8.2, margin = margin(b = 1)),
      plot.subtitle = element_text(hjust = 0.5, size = 7, margin = margin(b = 3)),
      axis.title = element_text(size = 7.2),
      axis.text = element_text(size = 6.6, color = "black"),
      axis.text.x = element_text(lineheight = 0.9),
      axis.line = element_line(linewidth = 0.38),
      axis.ticks = element_line(linewidth = 0.34),
      plot.margin = margin(5, 6, 5, 6)
    )

  if (nrow(bracket_df) > 0) {
    p <- p +
      geom_segment(data = bracket_df, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.34) +
      geom_segment(data = bracket_df, aes(x = x1, xend = x1, y = y_tip, yend = y), inherit.aes = FALSE, linewidth = 0.34) +
      geom_segment(data = bracket_df, aes(x = x2, xend = x2, y = y_tip, yend = y), inherit.aes = FALSE, linewidth = 0.34) +
      geom_text(data = bracket_df, aes(x = (x1 + x2) / 2, y = y + 0.035 * y_span, label = label), inherit.aes = FALSE, size = 2.05)
  }

  p
}

make_ancestry_plot <- function(target) {
  table_file <- file.path(
    table_dir,
    paste0(target$lead_snp, "_1000g_ancestry_genotype_proportions.tsv")
  )
  if (!file.exists(table_file)) {
    stop("Missing packaged 1000 Genomes summary: ", table_file)
  }

  df <- read_tsv(table_file, show_col_types = FALSE)
  ref <- unique(df$ref)[1]
  alt <- unique(df$alt)[1]
  genotype_levels <- c(
    paste0(ref, ref, " (0 ALT)"),
    paste0(ref, alt, " (1 ALT)"),
    paste0(alt, alt, " (2 ALT)")
  )

  df <- df %>%
    mutate(
      super_pop = factor(super_pop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS")),
      ancestry_label = factor(ancestry_label, levels = c("AFR", "AMR", "EAS", "EUR", "SAS")),
      genotype_label = factor(genotype_label, levels = genotype_levels)
    )

  labels <- df %>% distinct(ancestry_label, top_label)

  ggplot(df, aes(x = ancestry_label, y = proportion, fill = genotype_label)) +
    geom_col(width = 0.66, color = "#edf3f8", linewidth = 0.22) +
    geom_text(data = labels, aes(x = ancestry_label, y = 1.045, label = top_label), inherit.aes = FALSE, size = 1.65, lineheight = 0.86) +
    scale_x_discrete(expand = expansion(add = 0.55)) +
    scale_y_continuous(labels = number_format(accuracy = 0.1), breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(
      values = setNames(c("#c6dbef", "#6baed6", "#08519c"), genotype_levels),
      limits = genotype_levels,
      drop = FALSE,
      name = paste0(target$lead_snp, " REF=", ref, "/ALT=", alt),
      guide = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    coord_cartesian(ylim = c(0, 1.09), clip = "off") +
    labs(x = "Genetic ancestry", y = "Proportion of samples") +
    theme_classic(base_size = 7.4) +
    theme(
      axis.title = element_text(size = 7.1),
      axis.text = element_text(size = 6.4, color = "black"),
      axis.line = element_line(linewidth = 0.38),
      axis.ticks = element_line(linewidth = 0.34),
      legend.position = "bottom",
      legend.title = element_text(size = 5.8),
      legend.text = element_text(size = 5.7),
      legend.key.size = unit(0.18, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-3, 0, 0, 0),
      plot.margin = margin(10, 8, 7, 8)
    )
}

save_panel <- function(plot, file_stub) {
  ggsave(file.path(plot_dir, paste0(file_stub, ".png")), plot, width = 3, height = 3, units = "in", dpi = 600, bg = "white")
  ggsave(file.path(plot_dir, paste0(file_stub, ".pdf")), plot, width = 3, height = 3, units = "in", device = "pdf", bg = "white")
}

write_gtex_summary <- function(target) {
  summary <- load_gtex_target(target)$summary
  write_tsv(summary, file.path(table_dir, paste0(target$lead_snp, "_", target$target_gene, "_", target$gtex_tissue, "_gtex_summary.tsv")))
  invisible(summary)
}
