#!/usr/bin/env Rscript

source(file.path("scripts", "08_figS4_eqtl_ancestry", "R", "panel_helpers.R"))

# I show the thyroid eQTL pattern for the hypothyroidism FOXE1 lead SNP.
target <- target_row("rs10983761", "FOXE1", "Thyroid")
plot <- make_gtex_plot(target)
save_panel(plot, "rs10983761_FOXE1_Thyroid_gtex_eqtl")
write_gtex_summary(target)
