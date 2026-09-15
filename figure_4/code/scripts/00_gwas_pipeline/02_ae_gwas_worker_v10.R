#!/usr/bin/env Rscript
# =============================================================================
# AE GWAS v10 - WORKER (sensitivity: adjusting for IMPACT panel version)
# Identical to v9 worker except: reads v10 config, adds panel dummies to covariates
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(stringr); library(survival); library(SPACox)
})
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
cat("=============================================================================\n")
cat("   AE GWAS v10 - WORKER | Task:", task_id, "\n")
cat("   Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")
PROJECT_DIR <- "/data1/reznike/walserr_gwas"
CACHE_DIR   <- file.path(PROJECT_DIR, "cache_ae_gwas_v10")
config      <- readRDS(file.path(CACHE_DIR, "ae_gwas_config.rds"))
analysis_list <- config$analysis_list
if (task_id >= nrow(analysis_list)) {
  cat("Task ID", task_id, "exceeds number of analyses (", nrow(analysis_list), ")\n")
  quit(save = "no", status = 0)
}
row_idx     <- task_id + 1
grade_type  <- as.character(analysis_list$grade_type[row_idx])
event_def   <- as.character(analysis_list$event_def[row_idx])
cohort_type <- as.character(analysis_list$cohort_type[row_idx])
ae_name     <- as.character(analysis_list$ae_name[row_idx])
RESULTS_DIR <- as.character(analysis_list$results_dir[row_idx])
GENO_PREFIX <- config$geno_prefix
MIN_EVENTS  <- config$min_events
PC_COLS     <- config$pc_cols
GW_SIG      <- 5e-8
SUGGESTIVE  <- 1e-5
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
cat("  Analysis:", analysis_list$analysis_name[row_idx], "\n")
cat("  Grade:", grade_type, "| Event def:", event_def, "| Cohort:", cohort_type, "| AE:", ae_name, "\n")
cat("\n=== LOADING GENOTYPES ===\n")
fam <- fread(paste0(GENO_PREFIX, ".fam"), header = FALSE)
setnames(fam, c("FID", "IID", "PID", "MID", "SEX_PLINK", "PHENO"))
geno_ids <- fam$IID
cat("  Genotype samples:", format(nrow(fam), big.mark = ","), "\n")
cat("\n=== SELECTING COHORT ===\n")
cohort <- switch(cohort_type,
  ICI      = copy(config$cohort_ICI),
  NonICI   = copy(config$cohort_NonICI),
  Combined = copy(config$cohort_Combined)
)
cat("  Cohort size:", format(nrow(cohort), big.mark = ","), "\n")
cat("\n=== FINDING AE EVENTS ===\n")
pheno <- copy(cohort)
is_composite <- ae_name %in% config$ae_types_composite
llm <- config$llm_batch
if (is_composite) {
  if (ae_name == "thyroiditis") { components <- c("hyperthyroidism", "hypothyroidism")
  } else if (ae_name == "endocrine_irae") { components <- c("adrenal_insufficiency", "hyperthyroidism", "hypothyroidism")
  } else if (ae_name == "any_ae") { components <- config$ae_types_base }
  llm[, ae_positive := as.integer(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = components]
} else {
  llm[, ae_positive := as.integer(get(ae_name) > 0)]
}
llm_pos <- llm[ae_positive == 1L & dmp_id %in% pheno$dmp_id,
                .(dmp_id, window_start = window_start_days_from_dx,
                  window_end = window_end_days_from_dx)]
llm_pos <- merge(llm_pos, pheno[, .(dmp_id, t0, censor_date)], by = "dmp_id")
ae_dates <- llm_pos[window_start >= t0 & window_start <= censor_date,
                    .(ae_date = min(window_start)), by = dmp_id]
pheno[, ae_date := as.numeric(NA)]
pheno[ae_dates, ae_date := i.ae_date, on = "dmp_id"]
pheno[, event := as.integer(!is.na(ae_date))]
pheno[, time_to_event := fifelse(event == 1L, as.numeric(ae_date - t0), as.numeric(censor_date - t0))]
pheno <- pheno[time_to_event > 0]
n_events <- sum(pheno$event)
cat("  Final N:", format(nrow(pheno), big.mark = ","), "| Events:", n_events, "\n")
if (n_events < MIN_EVENTS) {
  cat("  SKIPPING: Too few events (", n_events, "<", MIN_EVENTS, ")\n")
  fwrite(data.table(analysis=analysis_list$analysis_name[row_idx], status="skipped_low_events",
                    n=nrow(pheno), events=n_events),
         file.path(RESULTS_DIR, paste0(ae_name, "_summary_v10.csv")))
  quit(save = "no", status = 0)
}
cat("\n=== BUILDING COVARIATES ===\n")
cancer_dummy_cols <- config$cancer_dummy_cols
panel_dummy_cols <- if (cohort_type == "ICI") config$panel_dummy_cols_ici else config$panel_dummy_cols_combined
if (cohort_type == "ICI") {
  covar_cols <- c("Sex", "Age", "LOT", "contains_chemo", "contains_ctla4",
                  "contains_pd1", "contains_biologic", "contains_targeted",
                  cancer_dummy_cols, panel_dummy_cols, PC_COLS)
} else {
  covar_cols <- c("Sex", "Age", "LOT", "contains_chemo", "contains_ctla4",
                  "contains_pd1", "contains_biologic", "contains_targeted",
                  cancer_dummy_cols, panel_dummy_cols, PC_COLS)
}
missing_cols <- setdiff(covar_cols, names(pheno))
if (length(missing_cols) > 0) {
  cat("  WARNING: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  covar_cols <- intersect(covar_cols, names(pheno))
}
for (col in covar_cols) pheno[, (col) := as.numeric(get(col))]
pheno <- pheno[complete.cases(pheno[, ..covar_cols])]
cat("  After complete cases:", format(nrow(pheno), big.mark = ","), "\n")
n_events <- sum(pheno$event)
if (n_events < MIN_EVENTS) {
  cat("  SKIPPING: Too few events after covariate filter (", n_events, ")\n")
  fwrite(data.table(analysis=analysis_list$analysis_name[row_idx], status="skipped_low_events_post_filter",
                    n=nrow(pheno), events=n_events),
         file.path(RESULTS_DIR, paste0(ae_name, "_summary_v10.csv")))
  quit(save = "no", status = 0)
}
cat("\n=== ORDERING BY FAM FILE ===\n")
pheno <- merge(fam[, .(IID, fam_order = .I)], pheno, by = "IID")
setorder(pheno, fam_order)
pheno[, fam_order := NULL]
cat("  Phenotype rows after fam merge:", format(nrow(pheno), big.mark = ","), "\n")
n_events <- sum(pheno$event)
cat("  Events:", n_events, "\n")
drop_cols <- c()
for (col in covar_cols) {
  if (length(unique(pheno[[col]])) < 2) {
    cat("  Dropping zero-variance:", col, "\n")
    drop_cols <- c(drop_cols, col)
  }
}
covar_cols <- setdiff(covar_cols, drop_cols)
cat("  Final covariates (", length(covar_cols), "):", paste(covar_cols, collapse = ", "), "\n")
cat("  Panel dummies included:", paste(intersect(panel_dummy_cols, covar_cols), collapse = ", "), "\n")
cat("\n=== FITTING NULL MODEL ===\n")
formula_str <- paste0("Surv(time_to_event, event) ~ ", paste(covar_cols, collapse = " + "))
cat("  Formula:", formula_str, "\n")
pheno_df <- as.data.frame(pheno)
formula_obj <- as.formula(gsub("time_to_event", "pheno_df$time_to_event",
                               gsub("event\\)", "pheno_df$event)", formula_str)))
obj_null <- tryCatch({
  SPACox_Null_Model(formula_obj, data = pheno_df, pIDs = pheno_df$IID, gIDs = geno_ids)
}, error = function(e) {
  cat("  ERROR fitting null model:", e$message, "\n")
  fwrite(data.table(analysis=analysis_list$analysis_name[row_idx], status="null_model_failed",
                    n=nrow(pheno), events=n_events, error=e$message),
         file.path(RESULTS_DIR, paste0(ae_name, "_summary_v10.csv")))
  quit(save = "no", status = 1)
})
cat("  Null model fitted successfully\n")
cat("\n=== RUNNING SPACox.plink ===\n")
output_file <- file.path(RESULTS_DIR, paste0(ae_name, "_SPACox_v10.txt"))
if (file.exists(output_file)) file.remove(output_file)
start_time <- Sys.time()
tryCatch({
  SPACox.plink(obj_null, GENO_PREFIX, output_file)
}, error = function(e) {
  cat("  ERROR in SPACox.plink:", e$message, "\n")
  fwrite(data.table(analysis=analysis_list$analysis_name[row_idx], status="spacox_failed",
                    n=nrow(pheno), events=n_events, error=e$message),
         file.path(RESULTS_DIR, paste0(ae_name, "_summary_v10.csv")))
  quit(save = "no", status = 1)
})
runtime <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
cat("  Completed in", runtime, "minutes\n")
cat("\n=== PROCESSING RESULTS ===\n")
gwas <- fread(output_file)
p_col <- names(gwas)[grep("p.value", names(gwas), ignore.case = TRUE)[1]]
if (is.na(p_col)) p_col <- names(gwas)[ncol(gwas)]
gwas[, P := as.numeric(get(p_col))]
gwas <- gwas[!is.na(P) & P > 0 & P <= 1]
snp_col <- names(gwas)[grep("marker|snp", names(gwas), ignore.case = TRUE)[1]]
if (!is.na(snp_col) && snp_col != "SNP") setnames(gwas, snp_col, "SNP")
bim <- fread(paste0(GENO_PREFIX, ".bim"), header = FALSE)
setnames(bim, c("CHR", "SNP", "CM", "BP", "A1", "A2"))
gwas <- merge(gwas, bim[, .(SNP, CHR, BP)], by = "SNP", all.x = TRUE)
gwas <- gwas[!is.na(CHR) & !is.na(BP)]
cat("  Valid SNPs:", format(nrow(gwas), big.mark = ","), "\n")
chisq <- qchisq(1 - gwas$P, 1)
lambda_gc <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
cat("  Lambda GC:", round(lambda_gc, 4), "\n")
n_gw <- sum(gwas$P < GW_SIG)
n_sugg <- sum(gwas$P < SUGGESTIVE)
cat("  Genome-wide significant (p<5e-8):", n_gw, "\n")
cat("  Suggestive (p<1e-5):", n_sugg, "\n")
top1000 <- gwas[order(P)][1:min(1000, nrow(gwas))]
fwrite(top1000, file.path(RESULTS_DIR, paste0(ae_name, "_top1000_v10.csv")))
if (n_sugg > 0) {
  fwrite(gwas[P < SUGGESTIVE][order(P)], file.path(RESULTS_DIR, paste0(ae_name, "_suggestive_v10.csv")))
}
system(paste("gzip -f", output_file))
summary_dt <- data.table(
  analysis=analysis_list$analysis_name[row_idx], grade_type=grade_type, event_def=event_def,
  cohort_type=cohort_type, ae_name=ae_name, n_samples=nrow(pheno), n_events=n_events,
  lambda_gc=round(lambda_gc,4), n_gw_sig=n_gw, n_suggestive=n_sugg, runtime_min=runtime,
  status="completed", covariates=paste(covar_cols, collapse=";"))
fwrite(summary_dt, file.path(RESULTS_DIR, paste0(ae_name, "_summary_v10.csv")))
cat("\n=== TASK COMPLETE ===\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
