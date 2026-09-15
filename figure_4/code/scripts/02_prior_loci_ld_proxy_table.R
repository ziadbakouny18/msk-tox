# Prior-loci (Groha et al. 2022) LD-proxy supplementary table for Figure 4.
#
# The manuscript tests the three published pooled-irAE loci against MSK-Tox and
# reports a null in aggregate plus two nominal, toxicity-specific hits. No source
# table exists. This builds it from the assembled Gusev lookup
# (gusev_nature_medicine_validation/msk_validation/results/).
#
# The null IS the headline here, so every (locus x AE x grade x cohort) cell gets
# a row. A table of just the two hits would invert the argument.
#
# Two things worth knowing before reading the output:
#   1. The two nominal hits sit at grade 2+, and the r2 that reproduces the
#      manuscript's rs149813236 is 0.75 -- i.e. below the r2 >= 0.8 the Methods
#      claim. The threshold used here is 0.6 (the panel these proxies came from);
#      the Methods sentence needs fixing, not the table. r2 is 1000G EUR.
#   2. SPACox returns a score test, not a Cox fit, so there is no HR to report.
#      Stat/Var degrades badly for rare large-effect variants (rs12705907 implies
#      HR 7.12 against 1.95 from the actual Cox fit), so this table reports z and
#      direction and deliberately does NOT report an HR.
#
# Writes: results/supplementary_table_prior_loci_ld_proxy.csv (collapsed)
#         results/prior_loci_all_candidate_proxies.csv         (full audit trail)

library(tidyverse)

root <- "/Users/elbakoz/Library/CloudStorage/Box-Box/Ziad El Bakouny (Collaborate)/2025/7_ae_gwas/gwas_transfer_package/scripts/ziad"
out  <- file.path(root, "manuscript_supplementary_tables", "results")

R2_MIN    <- 0.6     # what the quoted proxies actually satisfy (see note 1)
P_NOMINAL <- 0.05

# ---------------------------------------------------------------------------
# 0. The three published loci
# ---------------------------------------------------------------------------
# Groha S, et al. Nat Med 28:2584-2591 (2022), ref 43. HR is per the ALT allele
# of the paper's "Ref>Alt" and all three are risk-increasing in discovery, so the
# expected direction in MSK-Tox is positive on the ALT allele.
#
# The manuscript names only rs16906115 and rs113861051; rs75824728 (IL22RA1) is
# the third locus and must appear or the "three loci" sentence is unsupported.
groha <- tribble(
  ~prior_lead_rsid, ~nearest_gene, ~locus,    ~prior_chr, ~prior_pos, ~prior_ref, ~prior_alt,
  ~groha_profile_hr, ~groha_profile_p, ~groha_mgh_hr, ~groha_mgh_p, ~groha_trial_hr, ~groha_trial_p,
  "rs16906115",  "IL7",     "8q21.13", 8, 79712998, "G", "A", 2.0, 3.8e-9, 2.5, 1.9e-3, 1.2, 0.05,
  "rs75824728",  "IL22RA1", "1p36.11", 1, 24464626, "G", "A", 1.9, 8.4e-9, 0.90, 0.81, 0.77, 0.02,
  "rs113861051", "4p15",    "4p15.1",  4, 34448039, "A", "G", 2.0, 1.1e-8, 1.4, 0.43, 1.18, 0.14
) %>%
  mutate(groha_expected_dir = sign(log(groha_profile_hr)))

# ---------------------------------------------------------------------------
# 1. Read the assembled MSK-Tox lookup
# ---------------------------------------------------------------------------
raw <- read_csv(
  file.path(root, "gusev_nature_medicine_validation/msk_validation/results/msk_validation_gusev_snps.csv"),
  col_types = cols(.default = col_character())
) %>%
  mutate(across(c(SNP_POS, N_SAMPLES, N_EVENTS, MAC, MAF, INFO_SCORE,
                  Stat, Var, z, P_SPACOX, P_NORM, LAMBDA_GC,
                  r2_1000G_EUR, r2_1000G_EUR_AFR),
                ~ suppressWarnings(as.numeric(.x))))

dat <- raw %>%
  inner_join(groha, by = c("target_SNP" = "prior_lead_rsid")) %>%
  rename(prior_lead_rsid = target_SNP) %>%
  mutate(
    r2_to_lead   = if_else(target_role == "focal", 1, r2_1000G_EUR),
    well_imputed = INFO_FLAG != "poorly_imputed",
    # SPACox z is keyed to COUNTED_ALLELE; put every direction on the ALT allele
    # so it is comparable to Groha's.
    msk_dir_alt = case_when(
      is.na(z)                    ~ NA_real_,
      COUNTED_ALLELE_IS == "ALT"  ~ sign(z),
      COUNTED_ALLELE_IS == "REF"  ~ -sign(z),
      TRUE                        ~ NA_real_),
    direction_matches_groha = case_when(
      is.na(msk_dir_alt) ~ NA,
      msk_dir_alt == groha_expected_dir ~ TRUE,
      TRUE ~ FALSE),
    grade_threshold = recode(MSK_GRADE_TYPE,
                             all_grade  = "grade 1+ (all grade)",
                             grade2plus = "grade 2+",
                             grade3plus = "grade 3+"),
    cohort = recode(COHORT, ICI = "ICI-treated", NonICI = "non-ICI-treated",
                    Combined = "combined"))

# ---------------------------------------------------------------------------
# 2. Full candidate-proxy audit trail
# ---------------------------------------------------------------------------
audit <- dat %>%
  transmute(prior_lead_rsid, nearest_gene, cohort, grade_threshold, msk_ae = AE,
            role = target_role, variant_rsid = RSID, variant_id = SNP_ID,
            variant_chr = SNP_CHR, variant_pos = SNP_POS,
            msk_ref = MSK_REF, msk_alt = MSK_ALT,
            counted_allele = COUNTED_ALLELE, counted_allele_is = COUNTED_ALLELE_IS,
            r2_to_lead, r2_panel = "1000G EUR",
            maf = MAF, info_score = INFO_SCORE, info_flag = INFO_FLAG,
            n_samples = N_SAMPLES, n_events = N_EVENTS,
            z, p_spacox = P_SPACOX, msk_dir_alt, direction_matches_groha) %>%
  arrange(prior_lead_rsid, cohort, grade_threshold, msk_ae, desc(role), p_spacox)

write_csv(audit, file.path(out, "prior_loci_all_candidate_proxies.csv"))

# ---------------------------------------------------------------------------
# 3. Collapse to one row per locus x AE x grade x cohort
# ---------------------------------------------------------------------------
lead_rows <- dat %>%
  filter(target_role == "focal") %>%
  group_by(prior_lead_rsid, COHORT, MSK_GRADE_TYPE, AE) %>%
  slice_min(P_SPACOX, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(prior_lead_rsid, COHORT, MSK_GRADE_TYPE, AE,
            lead_tested_in_msk = TRUE,
            lead_maf = MAF, lead_info = INFO_SCORE,
            lead_z = z, lead_p = P_SPACOX,
            lead_dir_alt = msk_dir_alt,
            lead_dir_matches_groha = direction_matches_groha)

proxy_pool <- dat %>%
  filter(target_role == "proxy", well_imputed,
         !is.na(r2_to_lead), r2_to_lead >= R2_MIN)

best_by_p <- proxy_pool %>%
  group_by(prior_lead_rsid, COHORT, MSK_GRADE_TYPE, AE) %>%
  slice_min(P_SPACOX, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(prior_lead_rsid, COHORT, MSK_GRADE_TYPE, AE,
            best_proxy_rsid = RSID, best_proxy_id = SNP_ID,
            best_proxy_r2 = r2_to_lead,
            best_proxy_maf = MAF, best_proxy_info = INFO_SCORE,
            best_proxy_counted_allele = COUNTED_ALLELE,
            best_proxy_z = z, best_proxy_p = P_SPACOX,
            best_proxy_dir_alt = msk_dir_alt,
            best_proxy_dir_matches_groha = direction_matches_groha,
            n_samples = N_SAMPLES, n_events = N_EVENTS)

# The manuscript quotes rs149813236 for grade 2+ pneumonitis (p = 2.4e-3), which
# is INFO 0.271 -- below the 0.3 cut, i.e. the same INFO trap that bit the DFCI
# proxies. Keep the unfiltered pick as its own columns so that number stays
# traceable, with its INFO visible next to it, and so the better-imputed
# alternative can be compared directly.
best_by_p_any_info <- dat %>%
  filter(target_role == "proxy", !is.na(r2_to_lead), r2_to_lead >= R2_MIN) %>%
  group_by(prior_lead_rsid, COHORT, MSK_GRADE_TYPE, AE) %>%
  slice_min(P_SPACOX, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(prior_lead_rsid, COHORT, MSK_GRADE_TYPE, AE,
            unfiltered_best_proxy_rsid = RSID,
            unfiltered_best_proxy_r2 = r2_to_lead,
            unfiltered_best_proxy_info = INFO_SCORE,
            unfiltered_best_proxy_info_flag = INFO_FLAG,
            unfiltered_best_proxy_p = P_SPACOX,
            unfiltered_best_proxy_dir_matches_groha = direction_matches_groha)

proxy_counts <- dat %>%
  filter(target_role == "proxy") %>%
  group_by(prior_lead_rsid, COHORT, MSK_GRADE_TYPE, AE) %>%
  summarise(n_proxies_available = n(),
            n_proxies_tested = sum(!is.na(r2_to_lead) & r2_to_lead >= R2_MIN & well_imputed),
            n_proxies_dropped_low_info =
              sum(!is.na(r2_to_lead) & r2_to_lead >= R2_MIN & !well_imputed),
            n_proxies_nominal =
              sum(!is.na(r2_to_lead) & r2_to_lead >= R2_MIN & well_imputed &
                    !is.na(P_SPACOX) & P_SPACOX < P_NOMINAL),
            n_proxies_dir_concordant =
              sum(!is.na(r2_to_lead) & r2_to_lead >= R2_MIN & well_imputed &
                    !is.na(direction_matches_groha) & direction_matches_groha),
            .groups = "drop")

tbl <- dat %>%
  distinct(prior_lead_rsid, nearest_gene, locus, prior_chr, prior_pos,
           prior_ref, prior_alt, groha_profile_hr, groha_profile_p,
           groha_mgh_hr, groha_mgh_p, groha_trial_hr, groha_trial_p,
           COHORT, MSK_GRADE_TYPE, AE, cohort, grade_threshold) %>%
  left_join(lead_rows,    by = c("prior_lead_rsid", "COHORT", "MSK_GRADE_TYPE", "AE")) %>%
  left_join(best_by_p,    by = c("prior_lead_rsid", "COHORT", "MSK_GRADE_TYPE", "AE")) %>%
  left_join(best_by_p_any_info,
                          by = c("prior_lead_rsid", "COHORT", "MSK_GRADE_TYPE", "AE")) %>%
  left_join(proxy_counts, by = c("prior_lead_rsid", "COHORT", "MSK_GRADE_TYPE", "AE")) %>%
  mutate(lead_tested_in_msk = coalesce(lead_tested_in_msk, FALSE),
         # The strongest evidence available in this cell, lead or proxy.
         reported_variant = case_when(
           !is.na(lead_p) & (is.na(best_proxy_p) | lead_p <= best_proxy_p) ~ "lead",
           !is.na(best_proxy_p) ~ "best proxy",
           TRUE ~ NA_character_),
         reported_rsid = if_else(reported_variant == "lead",
                                 prior_lead_rsid, best_proxy_rsid),
         reported_r2   = if_else(reported_variant == "lead", 1, best_proxy_r2),
         reported_z    = if_else(reported_variant == "lead", lead_z, best_proxy_z),
         reported_p    = if_else(reported_variant == "lead", lead_p, best_proxy_p),
         reported_dir_matches_groha = if_else(reported_variant == "lead",
                                              lead_dir_matches_groha,
                                              best_proxy_dir_matches_groha),
         # any_ae is the direct analogue of Groha's pooled-irAE phenotype; this is
         # the cell the "no association" claim rests on.
         is_pooled_analogue = AE == "any_ae",
         verdict = case_when(
           is.na(reported_p) ~ "not testable",
           reported_p < 5e-8 & reported_dir_matches_groha ~ "genome-wide significant, concordant",
           reported_p < P_NOMINAL & reported_dir_matches_groha ~ "nominal, direction concordant",
           reported_p < P_NOMINAL & !reported_dir_matches_groha ~ "nominal, direction discordant",
           TRUE ~ "no association")) %>%
  transmute(prior_lead_rsid, nearest_gene, locus,
            prior_chr, prior_pos, prior_effect_allele = prior_alt,
            prior_other_allele = prior_ref,
            groha_profile_hr, groha_profile_p, groha_mgh_hr, groha_mgh_p,
            groha_trial_hr, groha_trial_p,
            msk_ae = AE, is_pooled_analogue, grade_threshold, cohort,
            n_samples, n_events,
            lead_tested_in_msk, lead_maf, lead_info, lead_z, lead_p,
            lead_dir_matches_groha,
            n_proxies_available, n_proxies_tested, n_proxies_dropped_low_info,
            n_proxies_nominal, n_proxies_dir_concordant,
            best_proxy_rsid, best_proxy_r2, best_proxy_maf, best_proxy_info,
            best_proxy_z, best_proxy_p, best_proxy_dir_matches_groha,
            unfiltered_best_proxy_rsid, unfiltered_best_proxy_r2,
            unfiltered_best_proxy_info, unfiltered_best_proxy_info_flag,
            unfiltered_best_proxy_p, unfiltered_best_proxy_dir_matches_groha,
            reported_variant, reported_rsid, reported_r2, reported_z, reported_p,
            reported_dir_matches_groha, verdict) %>%
  arrange(prior_lead_rsid, msk_ae, grade_threshold, cohort)

write_csv(tbl, file.path(out, "supplementary_table_prior_loci_ld_proxy.csv"))

# ---------------------------------------------------------------------------
# 4. Shape only -- no rows to the console
# ---------------------------------------------------------------------------
cat("collapsed table rows:", nrow(tbl), " cols:", ncol(tbl), "\n")
cat("audit trail rows:", nrow(audit), "\n\n")
cat("verdict:\n"); print(table(tbl$verdict))
cat("\nverdict x locus:\n"); print(table(tbl$prior_lead_rsid, tbl$verdict))
cat("\npooled-irAE analogue (any_ae) cells by verdict:\n")
print(table(tbl$verdict[tbl$is_pooled_analogue]))
# Does the table reproduce the two numbers the manuscript already asserts?
cat("\nthe two manuscript-quoted cells, filtered vs unfiltered pick:\n")
tbl %>%
  filter((prior_lead_rsid == "rs16906115"  & msk_ae == "pneumonitis") |
         (prior_lead_rsid == "rs113861051" & msk_ae == "hypothyroidism"),
         grade_threshold == "grade 2+", cohort == "ICI-treated") %>%
  select(prior_lead_rsid, msk_ae,
         best_proxy_rsid, best_proxy_r2, best_proxy_info, best_proxy_p,
         unfiltered_best_proxy_rsid, unfiltered_best_proxy_r2,
         unfiltered_best_proxy_info, unfiltered_best_proxy_p) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nnominal cells (locus / AE / grade / cohort):\n")
tbl %>%
  filter(str_starts(verdict, "nominal")) %>%
  select(prior_lead_rsid, msk_ae, grade_threshold, cohort,
         reported_rsid, reported_r2, reported_p, verdict) %>%
  as.data.frame() %>%
  print(row.names = FALSE)
