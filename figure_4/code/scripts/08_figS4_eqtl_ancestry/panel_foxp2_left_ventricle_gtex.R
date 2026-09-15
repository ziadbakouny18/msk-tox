#!/usr/bin/env Rscript

source(file.path("scripts", "08_figS4_eqtl_ancestry", "R", "panel_helpers.R"))

# I include left ventricle as the strongest GTEx FOXP2 signal for this SNP.
target <- target_row("rs12705907", "FOXP2", "Heart_Left_Ventricle")
plot <- make_gtex_plot(target)
save_panel(plot, "rs12705907_FOXP2_Heart_Left_Ventricle_gtex_eqtl")
write_gtex_summary(target)
