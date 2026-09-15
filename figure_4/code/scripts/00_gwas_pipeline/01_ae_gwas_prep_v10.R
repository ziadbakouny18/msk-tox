#!/usr/bin/env Rscript
# =============================================================================
# AE GWAS v10 - PREP (sensitivity: adjusting for IMPACT panel version)
# Only 3 analyses: ICI pneumonitis start, ICI adrenal_insufficiency start,
#                   Combined hypothyroidism start (all_grade)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(stringr)})
cat("=============================================================================\n")
cat("   AE GWAS v10 - PREP (IMPACT panel sensitivity)\n")
cat("   Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================================\n\n")
PROJECT_DIR  <- "/data1/reznike/walserr_gwas"
CACHE_DIR_V9 <- file.path(PROJECT_DIR, "cache_ae_gwas_v9")
CACHE_DIR    <- file.path(PROJECT_DIR, "cache_ae_gwas_v10")
RESULTS_BASE <- file.path(PROJECT_DIR, "results_ae_gwas_v10")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_BASE, showWarnings = FALSE, recursive = TRUE)
cat("=== LOADING v9 CONFIG ===\n")
v9 <- readRDS(file.path(CACHE_DIR_V9, "ae_gwas_config.rds"))
cat("  v9 config loaded\n")
# Extract IMPACT panel version from IID (e.g., P-0105681-T01-IM7 -> IM7)
extract_panel <- function(iid) sub(".*-", "", iid)
for (cname in c("cohort_ICI", "cohort_NonICI", "cohort_Combined")) {
  v9[[cname]][, impact_panel := extract_panel(IID)]
}
cat("\n=== IMPACT PANEL DISTRIBUTION (ICI) ===\n")
print(v9$cohort_ICI[, .N, by = impact_panel][order(-N)])
cat("\n=== IMPACT PANEL DISTRIBUTION (Combined) ===\n")
print(v9$cohort_Combined[, .N, by = impact_panel][order(-N)])
# Build panel dummies (drop smallest panel as reference)
all_panels <- sort(unique(c(v9$cohort_ICI$impact_panel, v9$cohort_Combined$impact_panel)))
cat("\n  Unique panels:", paste(all_panels, collapse = ", "), "\n")
for (cname in c("cohort_ICI", "cohort_NonICI", "cohort_Combined")) {
  panel_counts <- v9[[cname]][, .N, by = impact_panel][order(-N)]
  ref_panel <- panel_counts$impact_panel[nrow(panel_counts)]
  panels_to_dummy <- panel_counts$impact_panel[panel_counts$impact_panel != ref_panel]
  for (p in panels_to_dummy) {
    col_name <- paste0("panel_", p)
    v9[[cname]][, (col_name) := as.integer(impact_panel == p)]
  }
  if (cname == "cohort_ICI") {
    panel_dummy_cols_ici <- paste0("panel_", panels_to_dummy)
    cat("  ICI panel dummies (ref=", ref_panel, "):", paste(panels_to_dummy, collapse=", "), "\n")
  }
  if (cname == "cohort_Combined") {
    panel_dummy_cols_combined <- paste0("panel_", panels_to_dummy)
    cat("  Combined panel dummies (ref=", ref_panel, "):", paste(panels_to_dummy, collapse=", "), "\n")
  }
}
# Build analysis list (3 analyses only)
analyses <- list(
  data.table(task_id=0, analysis_name="all_grade_ICI_pneumonitis_start",
             grade_type="all_grade", event_def="start", cohort_type="ICI",
             ae_name="pneumonitis", results_dir=file.path(RESULTS_BASE, "all_grade/start/ICI")),
  data.table(task_id=1, analysis_name="all_grade_ICI_adrenal_insufficiency_start",
             grade_type="all_grade", event_def="start", cohort_type="ICI",
             ae_name="adrenal_insufficiency", results_dir=file.path(RESULTS_BASE, "all_grade/start/ICI")),
  data.table(task_id=2, analysis_name="all_grade_Combined_hypothyroidism_start",
             grade_type="all_grade", event_def="start", cohort_type="Combined",
             ae_name="hypothyroidism", results_dir=file.path(RESULTS_BASE, "all_grade/start/Combined"))
)
analysis_list <- rbindlist(analyses)
for (d in unique(analysis_list$results_dir)) dir.create(d, showWarnings=FALSE, recursive=TRUE)
cat("\n=== ANALYSIS LIST ===\n")
print(analysis_list[, .(task_id, analysis_name)])
config <- list(
  cohort_ICI      = v9$cohort_ICI,
  cohort_NonICI   = v9$cohort_NonICI,
  cohort_Combined = v9$cohort_Combined,
  llm_batch       = v9$llm_batch,
  grade_data      = v9$grade_data,
  analysis_list   = analysis_list,
  cancer_dummy_cols = v9$cancer_dummy_cols,
  panel_dummy_cols_ici = panel_dummy_cols_ici,
  panel_dummy_cols_combined = panel_dummy_cols_combined,
  top_cancers     = v9$top_cancers,
  geno_prefix     = v9$geno_prefix,
  pc_cols         = v9$pc_cols,
  ae_types_base   = v9$ae_types_base,
  ae_types_composite = v9$ae_types_composite,
  ae_types_all    = v9$ae_types_all,
  min_events      = v9$min_events,
  censor_buffer   = v9$censor_buffer,
  version         = "v10"
)
saveRDS(config, file.path(CACHE_DIR, "ae_gwas_config.rds"))
cat("\n=== PREP COMPLETE ===\n")
cat("  Config saved to:", file.path(CACHE_DIR, "ae_gwas_config.rds"), "\n")
cat("  Analyses:", nrow(analysis_list), "\n")
