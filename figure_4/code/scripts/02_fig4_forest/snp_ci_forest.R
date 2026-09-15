# =============================================================================
# SNP irAE Cumulative Incidence & By-Cohort Forest Plots (COMPETING-RISK)
# -----------------------------------------------------------------------------
# Flat, run-top-to-bottom script. The data pipeline is inline so you can step
# through it: highlight a block and run it, inspect the objects it creates
# (master, dosage_mrn, surv, fd), then continue.
#
# To debug ONE analysis interactively:
#   1) run the libraries + CONFIG + SNP registry + SETTINGS + UTILITIES + the
#      two plot builders (everything above the MAIN PIPELINE).
#   2) set:   si <- 1            # which SNP (1..4)
#             snp <- SNPS[[si]]
#      then run the "PREP MASTER" and "DOSAGE" blocks.
#   3) set:   coh <- "ICI"       # or "non-ICI" / "combined"
#      then run the "SURVIVAL" block and inspect `surv`.
#
# Methods:
#   CI curves   = Aalen-Johansen CIF (cmprsk::cuminc), death competing.
#   Forest HRs  = Fine-Gray sHR (cmprsk::crr, death competing) AND cause-specific
#                 Cox HR (survival::coxph, death censored), per allele-dosage unit.
# Timeline (per anchor line, from REAL dates -- not t_*):
#   t0     = lot_start
#   censor = min(lot_end + 180, next_lot_start - 1, dod, last_fu)   (keep > t0)
#   event  = 1 irAE if ae_date in [t0, censor]; else 2 death if dod <= censor; else 0
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(lubridate)
library(cmprsk)
library(survival)
library(ggplot2)
library(patchwork)

# =============================================================================
# CONFIG  --  data reads
# =============================================================================
llm_hypoth  <- read_csv("/Users/elbakoz/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/0_research/1_projects/10_irAE/data/mskcc/line_of_therapy/tox_specific_lot_tables/llm84k_hypothyroidism_grade0_12dec2025.csv")
llm_pneum   <- read_csv("/Users/elbakoz/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/0_research/1_projects/10_irAE/data/mskcc/line_of_therapy/tox_specific_lot_tables/llm84k_pneumonitis_grade0_12dec2025.csv")
llm_adrenal <- read_csv("/Users/elbakoz/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/0_research/1_projects/10_irAE/data/mskcc/line_of_therapy/tox_specific_lot_tables/llm84k_adrenal_insufficiency_grade0_12dec2025.csv")

adrenal_snp              <- read_tsv("/Users/elbakoz/Library/CloudStorage/Box-Box/Ziad El Bakouny (Collaborate)/2025/7_ae_gwas/gwas_transfer_package/results/tomin_5_13_2026/snp_patientlevel/6_32537711.bcf.tsv")
pneumonitis_snp          <- read_tsv("/Users/elbakoz/Library/CloudStorage/Box-Box/Ziad El Bakouny (Collaborate)/2025/7_ae_gwas/gwas_transfer_package/results/tomin_5_13_2026/snp_patientlevel/7_113403267.bcf.tsv")
hypothyroidism_esrp1_snp <- read_tsv("/Users/elbakoz/Library/CloudStorage/Box-Box/Ziad El Bakouny (Collaborate)/2025/7_ae_gwas/gwas_transfer_package/results/tomin_5_13_2026/snp_patientlevel/8_95636867.bcf.tsv")
hypothyroidism_foxe1_snp <- read_tsv("/Users/elbakoz/Library/CloudStorage/Box-Box/Ziad El Bakouny (Collaborate)/2025/7_ae_gwas/gwas_transfer_package/results/tomin_5_13_2026/snp_patientlevel/9_100553957.bcf.tsv")

# Each SNP: which dosage table + which phenotype LOT table it pairs with.
# forest_xlim / forest_breaks: hardcoded per SNP so the small forest axis looks clean
# (no crowding/overlap). Widen forest_xlim if a whisker ever clips.
SNPS <- list(
  list(rsid = "rs12705907",  gene = "",      snp = pneumonitis_snp,          master = "pneum",
       phenotype = "pneumonitis",           chr = 7, pos = 113403267, ref = "T", alt = "C", x_truncate_days = 500,
       forest_xlim = c(0.35, 7),    forest_breaks = c(0.5, 1, 2, 5)),
  list(rsid = "rs112463084", gene = "",      snp = adrenal_snp,              master = "adrenal",
       phenotype = "adrenal_insufficiency", chr = 6, pos = 32537711,  ref = "A", alt = "G", x_truncate_days = NA,
       forest_xlim = c(0.45, 1.4),  forest_breaks = c(0.5, 0.75, 1, 1.25)),
  list(rsid = "rs10983761",  gene = "FOXE1", snp = hypothyroidism_foxe1_snp, master = "hypoth",
       phenotype = "hypothyroidism",        chr = 9, pos = 100553957, ref = "A", alt = "C", x_truncate_days = NA,
       forest_xlim = c(0.65, 2.7),  forest_breaks = c(0.75, 1, 1.5, 2)),
  list(rsid = "rs112386546", gene = "ESRP1", snp = hypothyroidism_esrp1_snp, master = "hypoth",
       phenotype = "hypothyroidism",        chr = 8, pos = 95636867,  ref = "G", alt = "A", x_truncate_days = NA,
       forest_xlim = c(0.8, 2.7),   forest_breaks = c(1, 1.5, 2, 2.5))
)
MASTER_TABLES <- list(pneum = llm_pneum, adrenal = llm_adrenal, hypoth = llm_hypoth)

# =============================================================================
# SETTINGS
# =============================================================================
PHENO_DISPLAY <- c(pneumonitis = "Pneumonitis",
                   adrenal_insufficiency = "Adrenal Insufficiency",
                   hypothyroidism = "Hypothyroidism")
COHORTS       <- c("ICI", "non-ICI", "combined")
MIN_EVENTS    <- 5      # minimum irAE events to fit a model / draw a curve
CENSOR_BUFFER <- 180    # days added to lot_end before censoring (matches v9)
VERSION       <- "cr_v2"

# Aggregate RESULTS (plots + estimate tables; NO row-level data) live with the script.
OUT_DIR <- "/Users/elbakoz/Library/CloudStorage/Box-Box/Ziad El Bakouny (Collaborate)/2025/7_ae_gwas/gwas_transfer_package/scripts/ziad/snp_forest/snp_ci_forest_R_outputs"
# Row-level / patient-level data must NOT be written into the scripts folder.
# The per-patient survival tables stay in memory (surv_all). If you ever want to
# persist them, write to this intermediate data folder instead -- never OUT_DIR.
DATA_DIR <- "/Users/elbakoz/Library/CloudStorage/Box-Box/Ziad El Bakouny (Collaborate)/2025/7_ae_gwas/gwas_transfer_package/data/intermediate"

COLOR_LOW <- "#1f4e79"; COLOR_HIGH <- "#c0392b"          # CI-curve dosage groups
COLOR_DEATH_LOW <- "#7fa9cb"; COLOR_DEATH_HIGH <- "#e08a7a"
# Cohort colors for the forest (match the supplied legend; edit if not exact).
COLOR_ICI <- "#3D5A80"; COLOR_NONICI <- "#E15A46"; COLOR_COMBINED <- "#000000"

# =============================================================================
# UTILITIES  (small, vectorized -- the only helpers in the script)
# =============================================================================
standardize_mrn <- function(x) {
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    d <- str_extract(gsub("[\"'P\\-]|MSK", "", trimws(as.character(s))), "\\d+")
    if (is.na(d)) NA_character_ else sprintf("%08d", as.integer(d))
  }, character(1), USE.NAMES = FALSE)
}
to_date <- function(x) {                                   # dates are m/d/y e.g. 1/22/21
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  as.Date(suppressWarnings(parse_date_time(as.character(x),
    orders = c("m/d/y","m/d/Y","Y-m-d","Y-m-d H:M:S","m/d/Y H:M:S"), quiet = TRUE)))
}
as_flag    <- function(x) if (is.logical(x)) x else toupper(trimws(as.character(x))) %in% c("TRUE","1","Y","YES","T")
coerce_sex <- function(x) { if (is.numeric(x)) return(as.numeric(x))
  s <- toupper(trimws(as.character(x))); ifelse(s %in% c("1","M","MALE"), 1, ifelse(s %in% c("0","F","FEMALE"), 0, NA_real_)) }
id_tag     <- function(s) if (!is.null(s$rsid) && !is.na(s$rsid)) s$rsid else s$gene
snp_label  <- function(s) sprintf("%s (chr%d:%s %s>%s)", id_tag(s), s$chr,
                                  format(s$pos, big.mark = ",", scientific = FALSE), s$ref, s$alt)

# =============================================================================
# PLOT BUILDERS  (two functions -- you won't usually debug these line-by-line)
# =============================================================================
# Cumulative incidence curve for one cohort, split by dosage > cutoff.
build_ci_plot <- function(surv, snp, cutoff, coh) {
  grp     <- as.integer(surv$dosage > cutoff)
  pheno   <- PHENO_DISPLAY[[snp$phenotype]]
  x_trunc <- if (!is.na(snp$x_truncate_days)) snp$x_truncate_days else max(surv$time_days)
  ci <- cuminc(surv$time_days, surv$event, group = grp, cencode = 0)
  curve_df <- NULL; legend_lab <- list()
  for (gv in c(0, 1)) {
    sub <- surv[grp == gv, ]
    if (nrow(sub) < 5 || sum(sub$event == 1) < 1) next
    col   <- if (gv == 0) COLOR_LOW else COLOR_HIGH
    col_d <- if (gv == 0) COLOR_DEATH_LOW else COLOR_DEATH_HIGH
    base  <- if (gv == 0) sprintf("dosage <= %.1f", cutoff) else sprintf("dosage > %.1f", cutoff)
    lab   <- sprintf("%s  (N=%d, events=%d)", base, nrow(sub), sum(sub$event == 1))
    if (!is.null(ci[[paste(gv, 1)]])) {
      curve_df <- bind_rows(curve_df, tibble(time = ci[[paste(gv,1)]]$time, est = ci[[paste(gv,1)]]$est,
                                             series = lab, color = col, lty = "solid")); legend_lab[[lab]] <- col }
    if (!is.null(ci[[paste(gv, 2)]]))
      curve_df <- bind_rows(curve_df, tibble(time = ci[[paste(gv,2)]]$time, est = ci[[paste(gv,2)]]$est,
                                             series = paste0(base, " (death)"), color = col_d, lty = "dashed"))
  }
  if (is.null(curve_df)) return(NULL)
  curve_df <- filter(curve_df, time <= x_trunc)
  # Fine-Gray sHR for the carrier split (title annotation)
  hr <- NA; pval <- NA
  if (length(unique(grp)) == 2 && sum(surv$event == 1) >= MIN_EVENTS) {
    fgc <- tryCatch({ f <- crr(surv$time_days, surv$event, cov1 = matrix(grp, ncol = 1), failcode = 1, cencode = 0)
                      co <- summary(f)$coef; c(co[1, 2], co[1, ncol(co)]) }, error = function(e) c(NA, NA))
    hr <- fgc[1]; pval <- fgc[2]
  }
  hr_str <- if (!is.na(hr)) sprintf("sHR = %.2f", hr) else "sHR = NA"
  p_str  <- if (is.na(pval)) "p = NA" else if (pval < 1e-3) sprintf("p = %.2e", pval) else sprintf("p = %.3f", pval)
  p_ci <- ggplot(curve_df, aes(time, est, group = series, color = color, linetype = lty)) +
    geom_step(linewidth = 0.8) +
    scale_color_identity(guide = "legend", breaks = unname(unlist(legend_lab)), labels = names(legend_lab)) +
    scale_linetype_identity() +
    coord_cartesian(xlim = c(0, x_trunc), ylim = c(0, NA)) +
    labs(x = "Time from t = 0 (days)", y = sprintf("Cumulative incidence of %s", tolower(pheno)), color = NULL,
         title = sprintf("%s -- %s | %s cohort | cutoff > %.1f | %s %s (Fine-Gray)",
                         snp_label(snp), pheno, coh, cutoff, hr_str, p_str),
         caption = "Dashed = death (competing risk)") +
    theme_classic(base_size = 11) +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_blank(), plot.title = element_text(size = 10))
  ticks <- scales::breaks_pretty()(c(0, x_trunc)); ticks <- ticks[ticks >= 0 & ticks <= x_trunc]
  risk <- expand.grid(t = ticks, g = c(0, 1)) |> rowwise() |>
    mutate(n = sum(surv$time_days >= t & grp == g), col = if (g == 0) COLOR_LOW else COLOR_HIGH) |> ungroup()
  p_risk <- ggplot(risk, aes(t, factor(g, levels = c(1, 0)), label = n, color = col)) +
    geom_text(size = 3) + scale_color_identity() +
    scale_y_discrete(labels = function(b) ifelse(b == "0", sprintf("dosage <= %.1f", cutoff), sprintf("dosage > %.1f", cutoff))) +
    coord_cartesian(xlim = c(0, x_trunc)) + labs(x = NULL, y = NULL, title = "Number at risk") +
    theme_classic(base_size = 9) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
          plot.title = element_text(size = 9, face = "bold"))
  p_ci / p_risk + plot_layout(heights = c(5, 1))
}

# Compact forest (subdistribution HR only), 3 cohort rows colored by cohort.
# Sized for ~170 x 185 pt -- minimal ink, no numeric labels (those are in the CSV).
build_forest_plot <- function(snp, fd) {
  fd <- fd |> filter(method == "Fine-Gray") |>
    mutate(cohort = factor(cohort, levels = rev(COHORTS)))   # ICI at top
  xlim <- snp$forest_xlim       # hardcoded per SNP in the SNPS registry (covers its CI range)
  brks <- snp$forest_breaks     # hardcoded ticks chosen to look clean (no overlap)
  ggplot(fd, aes(HR, cohort, color = cohort)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, color = "grey55") +
    geom_errorbar(aes(xmin = LCI, xmax = UCI), orientation = "y", width = 0.22, linewidth = 0.5) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("ICI" = COLOR_ICI, "non-ICI" = COLOR_NONICI, "combined" = COLOR_COMBINED),
                       guide = "none") +
    scale_x_log10(breaks = brks) +
    coord_cartesian(xlim = xlim) +
    labs(x = "HR per allele dosage (95% CI)", y = NULL,
         title = sprintf("%s - %s", id_tag(snp), PHENO_DISPLAY[[snp$phenotype]])) +
    theme_classic(base_size = 6) +
    theme(plot.title   = element_text(size = 6.5, face = "bold"),
          axis.text    = element_text(size = 6, color = "black"),
          axis.title.x = element_text(size = 6),
          axis.line    = element_line(linewidth = 0.3),
          axis.ticks   = element_line(linewidth = 0.3),
          plot.margin  = margin(2, 4, 2, 2))
}

# #############################################################################
# MAIN PIPELINE  (flat; objects land in your environment as it runs)
# #############################################################################
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
surv_all <- list(); estimates_all <- list()

for (si in seq_along(SNPS)) {
  snp <- SNPS[[si]]
  cat(sprintf("\n===== [%d/%d] %s | %s =====\n", si, length(SNPS), snp_label(snp), snp$phenotype))

  ## ---- 1) PREP MASTER (parse dates; LOT number from "10+" etc.) ----
  mraw <- MASTER_TABLES[[snp$master]]; names(mraw) <- trimws(names(mraw))
  master <- tibble(
    MRN             = standardize_mrn(mraw$mrn),
    SAMPLE_ID       = trimws(as.character(mraw$impact_sample_id)),
    lot_num         = suppressWarnings(as.integer(str_extract(as.character(mraw$lot), "\\d+"))),
    contains_immuno = as_flag(mraw$contains_immuno),
    lot_start       = to_date(mraw$lot_start),
    lot_end         = to_date(mraw$lot_end),
    next_lot_start  = to_date(mraw$next_lot_start),
    dod             = to_date(mraw$dod),
    last_fu         = to_date(mraw$last_fu),
    ae_date         = to_date(mraw$ae_date),
    age             = suppressWarnings(as.numeric(mraw$age_at_lot_start)),
    sex             = coerce_sex(mraw$sex),
    cancer_type     = as.character(mraw$cancer_type)
  ) |> filter(!is.na(MRN), !is.na(lot_num))

  ## ---- 2) DOSAGE -> one value per MRN (via impact_sample_id) ----
  dose <- snp$snp; names(dose) <- trimws(names(dose)); up <- toupper(names(dose))
  id_col  <- names(dose)[which(up %in% c("SAMPLE_ID","IID","FID"))[1]]; if (is.na(id_col)) id_col <- names(dose)[1]
  dos_col <- tail(setdiff(names(dose), id_col), 1)
  cat("  dosage file cols -> id:", id_col, "| dosage:", dos_col, "\n")
  dosage_raw <- tibble(SAMPLE_ID = trimws(as.character(dose[[id_col]])),
                       dosage    = suppressWarnings(as.numeric(dose[[dos_col]]))) |> filter(!is.na(dosage))
  sample_map <- master |> distinct(MRN, SAMPLE_ID) |> filter(!is.na(SAMPLE_ID), SAMPLE_ID != "", SAMPLE_ID != "NA")
  dosage_mrn <- dosage_raw |> inner_join(sample_map, by = "SAMPLE_ID") |>
    group_by(MRN) |> summarise(dosage = mean(dosage), .groups = "drop")
  d <- dosage_mrn$dosage; aaf <- mean(d) / 2
  cat(sprintf("  N=%d MRNs | AAF(alt=%s)=%.4f | MAF=%.4f | dosage>1.0: %.1f%% | >1.5: %.1f%%\n",
              length(d), snp$alt, aaf, min(aaf, 1 - aaf), 100 * mean(d > 1.0), 100 * mean(d > 1.5)))

  ## ---- patient-level earliest phenotype date (used by every cohort) ----
  ae_map <- master |> filter(MRN %in% dosage_mrn$MRN, !is.na(ae_date)) |>
    group_by(MRN) |> summarise(ae_date = min(ae_date), .groups = "drop")

  forest_rows <- list()
  for (coh in COHORTS) {

    ## ---- 3) ANCHOR LINE + COMPETING-RISK SURVIVAL (real dates) ----
    pts <- master |> filter(MRN %in% dosage_mrn$MRN)
    if (coh == "ICI")     pts <- filter(pts, contains_immuno)
    if (coh == "non-ICI") pts <- filter(pts, !contains_immuno)
    anchor <- pts |> group_by(MRN) |> slice_min(lot_num, n = 1, with_ties = FALSE) |> ungroup()

    surv <- anchor |>
      transmute(MRN, lot = lot_num, t0 = lot_start, lot_end, next_lot_start, dod, last_fu, age, sex, cancer_type) |>
      mutate(censor_date = pmin(lot_end + CENSOR_BUFFER, next_lot_start - 1, dod, last_fu, na.rm = TRUE)) |>
      filter(!is.na(t0), !is.na(censor_date), censor_date > t0) |>
      left_join(ae_map, by = "MRN") |>
      mutate(
        has_ae    = !is.na(ae_date) & ae_date >= t0 & ae_date <= censor_date,
        event     = if_else(has_ae, 1L, if_else(!is.na(dod) & dod <= censor_date, 2L, 0L)),
        time_days = if_else(has_ae, as.numeric(ae_date - t0), as.numeric(censor_date - t0))) |>
      filter(time_days > 0) |>
      inner_join(dosage_mrn, by = "MRN")
    surv_all[[paste(id_tag(snp), coh)]] <- surv      # kept for inspection

    cat(sprintf("  [%-8s] n=%d  AE=%d  death=%d  censored=%d\n",
                coh, nrow(surv), sum(surv$event == 1), sum(surv$event == 2), sum(surv$event == 0)))
    if (nrow(surv) < 20 || sum(surv$event == 1) < MIN_EVENTS) { cat("    too few events; skipping cohort\n"); next }

    ## ---- 4) CI CURVES (cutoffs 1.0 and 1.5) ----
    for (cut in c(1.0, 1.5)) {
      p <- tryCatch(build_ci_plot(surv, snp, cut, coh), error = function(e) { cat("    CI error:", conditionMessage(e), "\n"); NULL })
      if (!is.null(p)) {
        f <- file.path(OUT_DIR, sprintf("%s_%s_%s_ci_cut%.1f_%s.pdf", id_tag(snp), snp$phenotype, coh, cut, VERSION))
        ggsave(f, p, width = 11, height = 7.5, device = "pdf"); cat("    wrote", basename(f), "\n")
      }
    }

    ## ---- 5) PER-DOSAGE-UNIT MODELS (Fine-Gray + cause-specific Cox) ----
    n_ae <- sum(surv$event == 1)
    if (length(unique(surv$dosage)) < 3) { cat("    no dosage variation; skipping models\n"); next }
    cov <- data.frame(dosage = surv$dosage)                          # covariate of interest FIRST
    if (sum(!is.na(surv$age)) >= max(20, 0.5 * nrow(surv))) cov$age <- ifelse(is.na(surv$age), median(surv$age, na.rm = TRUE), surv$age)
    if (sum(!is.na(surv$sex)) >= max(20, 0.5 * nrow(surv))) cov$sex <- ifelse(is.na(surv$sex), median(surv$sex, na.rm = TRUE), surv$sex)
    cov$lot <- as.numeric(surv$lot)   # adjust for line of therapy (constant in combined -> auto-dropped below)
    adj <- setdiff(names(cov), "dosage"); adj <- adj[vapply(adj, function(c) sd(cov[[c]]) > 1e-8, logical(1))]
    cov <- cov[, c("dosage", adj), drop = FALSE]

    fg <- tryCatch({
      f  <- crr(surv$time_days, surv$event, cov1 = as.matrix(cov), failcode = 1, cencode = 0)
      co <- summary(f)$coef; ci <- summary(f)$conf.int
      data.frame(method = "Fine-Gray", HR = co["dosage", 2], LCI = ci["dosage", 3], UCI = ci["dosage", 4], p = co["dosage", ncol(co)])
    }, error = function(e) { cat("    Fine-Gray failed:", conditionMessage(e), "\n"); NULL })

    cx <- tryCatch({
      dat <- cbind(data.frame(time = surv$time_days, status = as.integer(surv$event == 1)), cov)
      f   <- coxph(as.formula(paste("Surv(time, status) ~", paste(names(cov), collapse = " + "))), data = dat)
      s   <- summary(f)
      data.frame(method = "Cox", HR = s$conf.int["dosage", "exp(coef)"], LCI = s$conf.int["dosage", "lower .95"],
                 UCI = s$conf.int["dosage", "upper .95"], p = s$coefficients["dosage", "Pr(>|z|)"])
    }, error = function(e) { cat("    Cox failed:", conditionMessage(e), "\n"); NULL })

    rows <- bind_rows(fg, cx)
    if (!is.null(rows) && nrow(rows)) { rows$cohort <- coh; rows$n <- nrow(surv); rows$events <- n_ae; forest_rows[[coh]] <- rows }
  }

  ## ---- 6) FOREST (3 cohorts; FG vs Cox) + estimates CSV ----
  if (length(forest_rows)) {
    fd <- bind_rows(forest_rows); estimates_all[[id_tag(snp)]] <- fd
    cat("\n  per-dosage-unit estimates:\n")
    print(as.data.frame(fd |> transmute(cohort, method, HR = round(HR, 3), LCI = round(LCI, 3),
                                        UCI = round(UCI, 3), p = signif(p, 3), n, events)), row.names = FALSE)
    ggsave(file.path(OUT_DIR, sprintf("%s_%s_forest_bycohort_%s.pdf", id_tag(snp), snp$phenotype, VERSION)),
           build_forest_plot(snp, fd), width = 170/72, height = 185/72, units = "in", device = "pdf")  # 170 x 185 pt
    write_csv(fd |> select(cohort, method, HR, LCI, UCI, p, n, events),
              file.path(OUT_DIR, sprintf("%s_%s_estimates_%s.csv", id_tag(snp), snp$phenotype, VERSION)))
    cat("  wrote forest PDF + estimates CSV\n")
  }
}

cat(sprintf("\nDone. Outputs in: %s\n", OUT_DIR))
# Inspect afterwards:  surv_all[["rs112463084 ICI"]]  |  estimates_all[["rs112463084"]]
