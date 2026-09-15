# DFCI replication supplementary table for the Figure 4 / GWAS section.
#
# The manuscript asserts specific DFCI replication HRs and P values but there is
# no source table anywhere on disk. This builds it from Zeyun's authoritative
# long file (msk_validation_long2.tsv, arrived 2026-08-26).
#
# One row per (MSK-Tox lead variant x DFCI phenotype x severity x ancestry).
# Failures get explicit rows -- as written, a reader only sees the two successes.
#
# Three things that are easy to get backwards, all handled below:
#   1. DFCI HR is per PROFILE_REF (COUNTED_ALLELE == PROFILE_REF in all 1,209
#      analysed rows). The manuscript reports HRs per ALT (Fig 4). So the
#      per-ALT column is the RECIPROCAL of the file's HR, with the CI inverted
#      and swapped. Asserted, not assumed.
#   2. Zeyun ignored INFO when assembling proxy sets, so picking the best proxy
#      by P lands on garbage (rs72751583, INFO 0.034, HR 7.98). Require
#      INFO >= 0.3 before ranking.
#   3. The 424 AFR rows have status NA and no fitted Cox model. They are neither
#      replication nor non-replication -- dropped, not counted as null.
#
# Writes: results/supplementary_table_dfci_gwas_replication.csv (collapsed)
#         results/dfci_replication_all_candidate_proxies.csv    (full audit trail)

library(tidyverse)

root <- "/Users/elbakoz/Library/CloudStorage/Box-Box/Ziad El Bakouny (Collaborate)/2025/7_ae_gwas/gwas_transfer_package/scripts/ziad"
out  <- file.path(root, "manuscript_supplementary_tables", "results")

R2_MIN   <- 0.8    # Methods: proxies at r2 >= 0.8 within +/- 500 kb
INFO_MIN <- 0.3    # keeps the INFO trap out of the "best proxy" pick
P_NOMINAL <- 0.05  # Methods: replication nominal at P < 0.05

# ---------------------------------------------------------------------------
# 0. The four MSK-Tox lead variants
# ---------------------------------------------------------------------------
# ref -> alt as used by the Fig 4 forest plots (per-ALT dosage). The long file
# identifies each locus only by target_CHR/target_POS, so this is the join key.
leads <- tribble(
  ~lead_rsid,      ~gene,      ~cytoband,  ~lead_chr, ~lead_pos,  ~msk_ref, ~msk_alt, ~discovery_ae,
  "rs10983761",    "FOXE1",    "9q22.33",  9,         100553957,  "A",      "C",      "Hypothyroidism",
  "rs112386546",   "ESRP1",    "8q22.1",   8,         95636867,   "G",      "A",      "Hypothyroidism",
  "rs12705907",    "FOXP2",    "7q31.1",   7,         113403267,  "T",      "C",      "Pneumonitis",
  "rs112463084",   "HLA-DRB1", "6p21.32",  6,         32537711,   "A",      "G",      "Adrenal insufficiency"
)

# ---------------------------------------------------------------------------
# 1. Read the DFCI long file and clean it up
# ---------------------------------------------------------------------------
long <- read_tsv(
  file.path(root, "dfci_zeyun_validation/8_15_2026_update/results/msk_validation_long2.tsv"),
  col_types = cols(.default = col_character())
) %>%
  mutate(across(c(SNP_POS, N_SAMPLES, N_GENO, N_EVENTS, MAC, MAF, INFO_SCORE,
                  HR, LOG_HR, SE_LOG_HR, P_COXPH, r2_PROFILE, r2_1000G),
                ~ suppressWarnings(as.numeric(.x))),
         target_CHR = as.integer(target_CHR),
         target_POS = as.integer(target_POS))

# The AFR set carries variant metadata but no fitted model. Drop it here and say
# so in the legend, rather than letting it look like 424 null results.
# NB read_tsv turns the literal "NA" status into a real NA, so test for that
# rather than for the string.
n_unfitted <- sum(is.na(long$status) | long$status != "OK")
long <- long %>% filter(!is.na(status), status == "OK")

# Guard the orientation assumption rather than trusting it.
stopifnot(all(long$COUNTED_ALLELE == long$PROFILE_REF))

long <- long %>%
  # HR_CI arrives as "0.87 (0.79-0.96)"; pull the bounds out.
  mutate(ci_lo_ref = as.numeric(str_match(HR_CI, "\\(([0-9.]+)-")[, 2]),
         ci_hi_ref = as.numeric(str_match(HR_CI, "-([0-9.]+)\\)")[, 2]),
         # Per ALT = reciprocal, bounds inverted and swapped.
         hr_alt    = 1 / HR,
         ci_lo_alt = 1 / ci_hi_ref,
         ci_hi_alt = 1 / ci_lo_ref,
         # r2 comes from PROFILE in-sample and/or 1000G; take the larger and
         # record which panel supplied it.
         r2_used   = pmax(r2_PROFILE, r2_1000G, na.rm = TRUE),
         r2_panel  = case_when(
           is.na(r2_PROFILE) & is.na(r2_1000G)                  ~ NA_character_,
           is.na(r2_1000G)                                      ~ "PROFILE",
           is.na(r2_PROFILE)                                    ~ "1000G",
           r2_PROFILE >= r2_1000G                               ~ "PROFILE",
           TRUE                                                 ~ "1000G"),
         well_imputed = !is.na(INFO_SCORE) & INFO_SCORE >= INFO_MIN) %>%
  inner_join(leads, by = c("target_CHR" = "lead_chr", "target_POS" = "lead_pos"))

# ---------------------------------------------------------------------------
# 2. Full candidate-proxy audit trail
# ---------------------------------------------------------------------------
# Every proxy actually tested, so a reviewer can see the pick wasn't cherry-run.
audit <- long %>%
  transmute(lead_rsid, gene, cytoband, ancestry = ANC, dfci_phenotype = AE,
            severity = SEVERITY, role = target_role,
            variant_rsid = RSID, variant_id = SNP_ID,
            variant_chr = SNP_CHR, variant_pos = SNP_POS,
            profile_ref = PROFILE_REF, profile_alt = PROFILE_ALT,
            counted_allele = COUNTED_ALLELE,
            maf = MAF, info_score = INFO_SCORE, well_imputed,
            r2_to_lead = r2_used, r2_panel, ld_source = LD_source,
            n_samples = N_SAMPLES, n_events = N_EVENTS,
            hr_per_ref = HR, hr_per_alt = hr_alt,
            ci_lo_per_alt = ci_lo_alt, ci_hi_per_alt = ci_hi_alt,
            p_value = P_COXPH) %>%
  arrange(lead_rsid, ancestry, dfci_phenotype, severity, desc(role), p_value)

write_csv(audit, file.path(out, "dfci_replication_all_candidate_proxies.csv"))

# ---------------------------------------------------------------------------
# 3. Collapse to one row per lead x phenotype x severity x ancestry
# ---------------------------------------------------------------------------
# The lead variant itself, where DFCI could test it at all.
lead_rows <- long %>%
  filter(target_role == "focal") %>%
  group_by(lead_rsid, ANC, AE, SEVERITY) %>%
  slice_min(P_COXPH, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(lead_rsid, ANC, AE, SEVERITY,
            lead_tested = TRUE,
            lead_info   = INFO_SCORE,
            lead_maf    = MAF,
            lead_hr_per_alt = hr_alt,
            lead_ci_lo  = ci_lo_alt,
            lead_ci_hi  = ci_hi_alt,
            lead_p      = P_COXPH,
            lead_well_imputed = well_imputed)

# Candidate proxies: r2 >= 0.8 and imputed well enough to believe.
proxy_pool <- long %>%
  filter(target_role == "proxy", well_imputed,
         !is.na(r2_used), r2_used >= R2_MIN)

# Methods define the reported proxy as the strongest DFCI association.
best_by_p <- proxy_pool %>%
  group_by(lead_rsid, ANC, AE, SEVERITY) %>%
  slice_min(P_COXPH, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(lead_rsid, ANC, AE, SEVERITY,
            best_proxy_rsid = RSID, best_proxy_id = SNP_ID,
            best_proxy_r2 = r2_used, best_proxy_r2_panel = r2_panel,
            best_proxy_info = INFO_SCORE, best_proxy_maf = MAF,
            best_proxy_counted_allele = COUNTED_ALLELE,
            best_proxy_alt_allele = PROFILE_ALT,
            best_proxy_hr_per_alt = hr_alt,
            best_proxy_ci_lo = ci_lo_alt, best_proxy_ci_hi = ci_hi_alt,
            best_proxy_p = P_COXPH,
            n_samples = N_SAMPLES, n_events = N_EVENTS)

# At HLA-DRB1 the r2-ranked and P-ranked picks differ by two orders of magnitude
# in P at near-identical r2, so report the r2-ranked one alongside it.
best_by_r2 <- proxy_pool %>%
  group_by(lead_rsid, ANC, AE, SEVERITY) %>%
  slice_max(r2_used, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(lead_rsid, ANC, AE, SEVERITY,
            highest_r2_proxy_rsid = RSID,
            highest_r2_proxy_r2 = r2_used,
            highest_r2_proxy_hr_per_alt = hr_alt,
            highest_r2_proxy_p = P_COXPH)

# How many proxies were on the table before ranking, and how many were thrown
# out for bad imputation -- both belong in the table, not a footnote.
proxy_counts <- long %>%
  filter(target_role == "proxy") %>%
  group_by(lead_rsid, ANC, AE, SEVERITY) %>%
  summarise(n_proxies_available = n(),
            n_proxies_r2_pass   = sum(!is.na(r2_used) & r2_used >= R2_MIN),
            n_proxies_tested    = sum(!is.na(r2_used) & r2_used >= R2_MIN & well_imputed),
            n_proxies_dropped_low_info =
              sum(!is.na(r2_used) & r2_used >= R2_MIN & !well_imputed),
            .groups = "drop")

tbl <- long %>%
  distinct(lead_rsid, gene, cytoband, discovery_ae, msk_ref, msk_alt, ANC, AE, SEVERITY) %>%
  left_join(lead_rows,    by = c("lead_rsid", "ANC", "AE", "SEVERITY")) %>%
  left_join(best_by_p,    by = c("lead_rsid", "ANC", "AE", "SEVERITY")) %>%
  left_join(best_by_r2,   by = c("lead_rsid", "ANC", "AE", "SEVERITY")) %>%
  left_join(proxy_counts, by = c("lead_rsid", "ANC", "AE", "SEVERITY")) %>%
  mutate(lead_tested = coalesce(lead_tested, FALSE),
         # n/events are per (phenotype x severity x ancestry) cell; take whichever
         # row carried them.
         n_samples = coalesce(n_samples, NA_real_),
         severity = recode(SEVERITY, mild = "grade 1+ (mild)",
                           mod = "grade 2+ (moderate)", sev = "grade 3+ (severe)"),

         # What could DFCI actually do with this locus?
         assessability = case_when(
           lead_tested & lead_well_imputed  ~ "lead directly assessable",
           lead_tested & !lead_well_imputed &
             !is.na(best_proxy_rsid)        ~ "lead poorly imputed; proxy-based",
           lead_tested & !lead_well_imputed ~ "lead poorly imputed; no usable proxy",
           !lead_tested & !is.na(best_proxy_rsid) ~ "lead absent; proxy-based",
           TRUE                             ~ "not assessable"),

         # The P value the manuscript should be quoting for this cell.
         reported_p  = if_else(assessability == "lead directly assessable",
                               lead_p, best_proxy_p),
         reported_hr = if_else(assessability == "lead directly assessable",
                               lead_hr_per_alt, best_proxy_hr_per_alt),
         discovery_dir = case_when(
           lead_rsid %in% c("rs10983761", "rs12705907") ~ 1,   # risk-increasing per ALT
           lead_rsid %in% c("rs112386546") ~ 1,
           lead_rsid == "rs112463084" ~ -1,                    # protective per ALT
           TRUE ~ NA_real_),
         dfci_dir = sign(log(reported_hr)),
         direction_concordant = case_when(
           is.na(dfci_dir) | is.na(discovery_dir) ~ NA,
           dfci_dir == discovery_dir ~ TRUE,
           TRUE ~ FALSE),

         replication_verdict = case_when(
           is.na(reported_p)                        ~ "not assessable",
           reported_p < 5e-8  & direction_concordant ~ "replicated (genome-wide)",
           reported_p < P_NOMINAL & direction_concordant ~ "replicated (nominal)",
           reported_p < P_NOMINAL & !direction_concordant ~ "nominal but direction discordant",
           TRUE                                     ~ "not replicated")) %>%
  transmute(lead_rsid, gene, cytoband, discovery_ae,
            msk_effect_allele = msk_alt, msk_other_allele = msk_ref,
            ancestry = ANC, dfci_phenotype = AE, severity,
            n_samples, n_events,
            assessability,
            lead_tested, lead_info, lead_maf,
            lead_hr_per_alt, lead_ci_lo, lead_ci_hi, lead_p,
            n_proxies_available, n_proxies_r2_pass, n_proxies_tested,
            n_proxies_dropped_low_info,
            best_proxy_rsid, best_proxy_r2, best_proxy_r2_panel,
            best_proxy_info, best_proxy_maf,
            best_proxy_effect_allele = best_proxy_alt_allele,
            best_proxy_hr_per_alt, best_proxy_ci_lo, best_proxy_ci_hi, best_proxy_p,
            highest_r2_proxy_rsid, highest_r2_proxy_r2,
            highest_r2_proxy_hr_per_alt, highest_r2_proxy_p,
            reported_hr_per_alt = reported_hr, reported_p,
            direction_concordant, replication_verdict) %>%
  arrange(lead_rsid, ancestry, dfci_phenotype, severity)

write_csv(tbl, file.path(out, "supplementary_table_dfci_gwas_replication.csv"))

# ---------------------------------------------------------------------------
# 4. Shape only -- no rows to the console
# ---------------------------------------------------------------------------
cat("AFR rows dropped (no fitted Cox model):", n_unfitted, "\n")
cat("collapsed table rows:", nrow(tbl), " cols:", ncol(tbl), "\n")
cat("audit trail rows:", nrow(audit), "\n\n")
cat("rows per lead:\n"); print(table(tbl$lead_rsid))
cat("\nassessability:\n"); print(table(tbl$assessability))
cat("\nreplication verdict:\n"); print(table(tbl$replication_verdict))
cat("\nverdict x lead:\n"); print(table(tbl$lead_rsid, tbl$replication_verdict))
