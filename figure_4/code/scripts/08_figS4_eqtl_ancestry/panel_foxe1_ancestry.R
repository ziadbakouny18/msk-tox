#!/usr/bin/env Rscript

source(file.path("scripts", "08_figS4_eqtl_ancestry", "R", "panel_helpers.R"))

# I show the 1000 Genomes genotype composition for the FOXE1 lead SNP.
target <- target_row("rs10983761", "FOXE1", "Thyroid")
plot <- make_ancestry_plot(target)
save_panel(plot, "rs10983761_FOXE1_1000g_ancestry")
