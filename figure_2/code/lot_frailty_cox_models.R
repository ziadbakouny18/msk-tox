#!/usr/bin/env Rscript
# ==============================================================================
# Line-of-therapy Cox models: shared-frailty (coxme, primary) + clustered-robust
# (coxph + cluster(), sensitivity) per toxicity.
#
# Input:  LOT_prevalence_S2E_cox_input.csv  (exported from the Python notebook)
#         One row per (patient, line of therapy), with:
#           mrn, lot, lot_group, age_at_lot_start, sex_clean, cancer_type_grp,
#           has_pd1_flag, has_ctla4_flag, contains_chemo, contains_hormone,
#           contains_biologic, contains_targeted,
#           {toxicity}_time_days, {toxicity}_event   for each of the 6 toxicities
#
# Adjustment covariates match the main figure panels (2C, 2D, 2E):
#   Treatment: has_pd1_flag, has_ctla4_flag, contains_chemo, contains_hormone,
#              contains_biologic, contains_targeted
#   Demographics: age_centered, sex_clean, cancer_type_grp
#
# Output: LOT_prevalence_S2E_frailty_cox_results.csv
#         One row per (toxicity, covariate, model_type), with HR, 95% CI, p-value,
#         plus frailty variance (theta) for the coxme rows.
# ==============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(coxme)
})

# ---- Paths (script lives in figure 2/scripts/, matching the Python notebook) ----
# Standard Rscript idiom for "this script's own directory" (works with `Rscript path/to/script.R`;
# if sourcing interactively instead, just set script_dir manually to the scripts/ folder).
# ---- Paths (assumes R's working directory is set to figure 2/scripts/, same
#      convention as the Python notebook's os.getcwd()) ----
# In RStudio: Session > Set Working Directory > To Source File Location (or setwd() manually).
# In a terminal: run `Rscript lot_frailty_cox_models.R` from inside the scripts/ folder.
script_dir  <- getwd()
results_dir <- file.path(dirname(script_dir), "results", "supp", "S2E_LOT_Prevalence")
input_path  <- file.path(results_dir, "LOT_prevalence_S2E_cox_input.csv")
output_path <- file.path(results_dir, "LOT_prevalence_S2E_frailty_cox_results.csv")

cat(sprintf("Working directory: %s\n", script_dir))
cat(sprintf("Looking for input at: %s\n", input_path))
if (!file.exists(input_path)) {
  stop(sprintf(
    "Input CSV not found at:\n  %s\nSet your R working directory to the 'figure 2/scripts' folder ",
    input_path
  ), call. = FALSE)
}
df <- read.csv(input_path, stringsAsFactors = FALSE)

TOXICITIES <- c("pneumonitis", "adrenal_insufficiency", "liver_toxicity",
                "colitis", "hyperthyroidism", "hypothyroidism")

# ---- Shared covariate prep (same for every toxicity) ----
# Matches the adjustment covariates used in main figure panels (2C, 2D, 2E)
df$lot_group  <- factor(df$lot_group, levels = c("1", "2", "3", "4+"))   # ref = LOT 1
df$sex_clean  <- factor(df$sex_clean, levels = c("Male", "Female"))       # ref = Male
df$cancer_type_grp <- factor(df$cancer_type_grp)
if ("Other" %in% levels(df$cancer_type_grp)) {
  df$cancer_type_grp <- relevel(df$cancer_type_grp, ref = "Other")       # ref = Other
}
df$age_centered <- df$age_at_lot_start - mean(df$age_at_lot_start, na.rm = TRUE)

# Treatment flags (ensure they're numeric 0/1)
df$has_pd1_flag <- as.integer(df$has_pd1_flag)
df$has_ctla4_flag <- as.integer(df$has_ctla4_flag)
df$contains_chemo <- as.integer(df$contains_chemo)
df$contains_hormone <- as.integer(df$contains_hormone)
df$contains_biologic <- as.integer(df$contains_biologic)
df$contains_targeted <- as.integer(df$contains_targeted)

fit_one_toxicity <- function(tox) {
  cat(sprintf("\n=== %s ===\n", tox))

  time_col  <- paste0(tox, "_time_days")
  event_col <- paste0(tox, "_event")
  d <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0, ]
  d$time  <- d[[time_col]]
  d$event <- d[[event_col]]

  # Model formula includes treatment flags for consistency with main figure panels
  form <- as.formula(
    "Surv(time, event) ~ lot_group + has_pd1_flag + has_ctla4_flag + contains_chemo + contains_hormone + contains_biologic + contains_targeted + age_centered + sex_clean + cancer_type_grp"
  )

  # ---- Primary: shared Gamma/log-normal frailty Cox model (coxme = Gaussian frailty
  #      on the log-hazard scale, the standard R equivalent for this design) ----
  frailty_form <- update(form, . ~ . + (1 | mrn))
  fit_frailty <- tryCatch(coxme(frailty_form, data = d), error = function(e) e)
  if (inherits(fit_frailty, "error")) {
    cat("  retrying with explicit starting variance (vinit)...\n")
    fit_frailty <- tryCatch(
      coxme(frailty_form, data = d, vinit = 0.5),
      error = function(e) e
    )
  }

  # ---- Sensitivity: clustered-robust-SE Cox (population-average HR, no variance
  #      component to estimate -- fallback if frailty is unstable for this toxicity) ----
  fit_clustered <- coxph(form, data = d, cluster = mrn)

  rows <- list()

  if (!inherits(fit_frailty, "error")) {
    cf   <- fixef(fit_frailty)
    vf   <- vcov(fit_frailty)
    se   <- sqrt(diag(vf))
    z    <- cf / se
    pval <- 2 * (1 - pnorm(abs(z)))
    # VarCorr() structure can vary slightly by coxme version; unlist+first element is robust
    theta <- as.numeric(unlist(VarCorr(fit_frailty)))[1]  # frailty variance on the log-hazard scale

    rows[["frailty"]] <- data.frame(
      toxicity = tox, model_type = "coxme_frailty",
      covariate = names(cf), beta = cf, se = se,
      HR = exp(cf), HR_ci_lower = exp(cf - 1.96 * se), HR_ci_upper = exp(cf + 1.96 * se),
      p_value = pval, frailty_theta = theta,
      n_obs = nrow(d), n_patients = length(unique(d$mrn)), n_events = sum(d$event),
      row.names = NULL
    )
  } else {
    cat("  coxme failed for this toxicity:", conditionMessage(fit_frailty), "\n")
  }

  cs   <- summary(fit_clustered)$coefficients
  rows[["clustered"]] <- data.frame(
    toxicity = tox, model_type = "clustered_robust_cox",
    covariate = rownames(cs), beta = cs[, "coef"], se = cs[, "robust se"],
    HR = exp(cs[, "coef"]),
    HR_ci_lower = exp(cs[, "coef"] - 1.96 * cs[, "robust se"]),
    HR_ci_upper = exp(cs[, "coef"] + 1.96 * cs[, "robust se"]),
    p_value = cs[, "Pr(>|z|)"], frailty_theta = NA,
    n_obs = nrow(d), n_patients = length(unique(d$mrn)), n_events = sum(d$event),
    row.names = NULL
  )

  do.call(rbind, rows)
}

all_results <- do.call(rbind, lapply(TOXICITIES, fit_one_toxicity))
write.csv(all_results, output_path, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", output_path))
