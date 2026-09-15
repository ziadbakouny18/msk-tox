# 01_scan_dfci_full_results.R
#
# One streaming pass over each DFCI (PROFILE) unfiltered SPA-Cox result file
# in results/full_gwas_results/, producing three combined tables:
#
#   qc_summary.tsv          per-trait N variants, lambda_GC, tail counts, min P.
#                           Zeyun re-ran everything after fixing a pipeline bug,
#                           so lambda_GC is the cheapest independent check that
#                           the re-run is well calibrated.
#   msk_lead_windows.tsv.gz every DFCI variant within +/-500 kb of each of our
#                           four reported MSK-Tox lead SNPs. This is the raw
#                           material for rebuilding the replication lookup that
#                           msk_validation_long.tsv was supposed to contain but
#                           which is not in the shared folder.
#   dfci_sig_snp_lookup.tsv the 15 DFCI genome-wide significant SNPs
#                           (dfci_sig_snps.tsv) looked up in every DFCI trait,
#                           for cross-phenotype context.
#
# Runtime is dominated by I/O: 16 files, ~5.6M variants each. Expect this to
# take a while; it is the only slow script in the folder.
#
# Everything read and written here is variant-level association summary
# statistics. No individual-level data is involved.

# ---- Packages ----
library(tidyverse)

# ---- Paths ----
# Run with the working directory set to the `analysis` folder, i.e. the parent
# of scripts/. In RStudio: Session > Set Working Directory > To Source File
# Location, then setwd("..").
analysis_dir <- normalizePath(".", mustWork = TRUE)
results_dir <- file.path(dirname(analysis_dir), "results")
full_dir <- file.path(results_dir, "full_gwas_results")
out_dir <- file.path(analysis_dir, "output")

stopifnot(dir.exists(file.path(analysis_dir, "scripts")), dir.exists(full_dir))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

window_bp <- 500000  # +/- 500 kb, matching the proxy window in results.docx

# ---- Reference inputs ----

# Our four MSK-Tox lead SNPs (GRCh37/hg19), from
# dfci_zeyun_validation/v10/collab_validation_summary.tsv.
msk_leads <- tribble(
  ~msk_lead_rsid,  ~lead_chr, ~lead_pos,  ~msk_lead_locus,
  "rs112463084",   6,         32537711,   "6p21.32 HLA-DRB1",
  "rs12705907",    7,         113403267,  "7q31.1 FOXP2",
  "rs112386546",   8,         95636867,   "8q22.1 ESRP1",
  "rs10983761",    9,         100553957,  "9q22.33 FOXE1"
)
glimpse(msk_leads)

# DFCI's own 15 genome-wide significant SNPs, which Zeyun asked us to look up.
dfci_sig <- read_tsv(file.path(results_dir, "dfci_sig_snps.tsv"),
                     show_col_types = FALSE)
stopifnot(all(c("RSID", "SNP_ID") %in% names(dfci_sig)))
glimpse(dfci_sig)

# The 16 traits, taken from the unfiltered file names.
traits <- list.files(full_dir, pattern = "_results\\.tsv\\.gz$") %>%
  discard(~ str_detect(.x, "_filtered_results")) %>%
  str_remove("^SPACox_") %>%
  str_remove("_results\\.tsv\\.gz$") %>%
  sort()
print(traits)
stopifnot(length(traits) == 16)

# Only the columns we need, so a 5.6M-row file stays manageable in memory.
scan_cols <- cols_only(
  SNP_ID = col_character(), MAF = col_double(), p.value.spa = col_double(),
  Stat = col_double(), Var = col_double(), z = col_double(),
  CHR = col_double(), BP = col_double(), gt_REF = col_character(),
  gt_ALT = col_character(), N = col_double(), RSID = col_character()
)

# ---- Scan ----
# One pass per trait. Accumulate into lists, bind once at the end.
qc_list <- vector("list", length(traits))
window_list <- vector("list", length(traits))
sig_list <- vector("list", length(traits))

for (i in seq_along(traits)) {
  trait <- traits[i]
  message(sprintf("[%2d/%d] %s", i, length(traits), trait))

  raw <- read_tsv(file.path(full_dir, sprintf("SPACox_%s_results.tsv.gz", trait)),
                  col_types = scan_cols, progress = FALSE)

  # Per-trait QC. lambda_GC uses the median p-value, converted back to a
  # 1-df chi-square, divided by the null median chi-square.
  median_p <- median(raw$p.value.spa, na.rm = TRUE)
  top_hit <- raw %>% slice_min(p.value.spa, n = 1, with_ties = FALSE)

  qc_list[[i]] <- tibble(
    trait = trait,
    ancestry = str_split_i(trait, "_", 1),
    severity = str_split_i(trait, "_", 2),
    adverse_event = str_remove(trait, "^[A-Z]+_[a-z]+_"),
    n_variants = nrow(raw),
    n_p_missing = sum(is.na(raw$p.value.spa)),
    n_samples = top_hit$N,
    median_p = median_p,
    lambda_gc = qchisq(median_p, df = 1, lower.tail = FALSE) / qchisq(0.5, df = 1),
    n_p_lt_5e3 = sum(raw$p.value.spa < 5e-3, na.rm = TRUE),
    n_p_lt_1e6 = sum(raw$p.value.spa < 1e-6, na.rm = TRUE),
    n_p_lt_5e8 = sum(raw$p.value.spa < 5e-8, na.rm = TRUE),
    min_p = top_hit$p.value.spa,
    min_p_rsid = top_hit$RSID,
    min_p_snp = top_hit$SNP_ID,
    min_p_maf = top_hit$MAF
  )

  # +/-500 kb windows around each MSK lead. The join keeps only the four
  # relevant chromosomes, then the filter trims to the window.
  window_list[[i]] <- raw %>%
    inner_join(msk_leads, by = c("CHR" = "lead_chr"), relationship = "many-to-many") %>%
    filter(abs(BP - lead_pos) <= window_bp) %>%
    mutate(trait = trait, dist_bp = BP - lead_pos) %>%
    select(trait, msk_lead_rsid, msk_lead_locus, lead_pos, dist_bp,
           SNP_ID, RSID, CHR, BP, gt_REF, gt_ALT, MAF, p.value.spa, z, N)

  # DFCI's own significant SNPs, in this trait.
  sig_list[[i]] <- raw %>%
    filter(SNP_ID %in% dfci_sig$SNP_ID) %>%
    mutate(trait = trait) %>%
    select(trait, SNP_ID, RSID, CHR, BP, gt_REF, gt_ALT, MAF, p.value.spa, z, N)

  rm(raw)
  gc(verbose = FALSE)
}

qc_summary <- bind_rows(qc_list) %>% arrange(ancestry, severity, adverse_event)
msk_lead_windows <- bind_rows(window_list)
dfci_sig_snp_lookup <- bind_rows(sig_list)

# ---- Inspect ----
print(qc_summary, n = Inf, width = Inf)
glimpse(msk_lead_windows)
glimpse(dfci_sig_snp_lookup)

# ---- Export ----
write_tsv(qc_summary, file.path(out_dir, "qc_summary.tsv"))
write_tsv(msk_lead_windows, file.path(out_dir, "msk_lead_windows.tsv.gz"))
write_tsv(dfci_sig_snp_lookup, file.path(out_dir, "dfci_sig_snp_lookup.tsv"))

message("wrote qc_summary.tsv, msk_lead_windows.tsv.gz, dfci_sig_snp_lookup.tsv to ", out_dir)
