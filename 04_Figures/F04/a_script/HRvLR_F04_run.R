#!/usr/bin/env Rscript
# HRvLR_F04_run.R — Run all F04 panels
# Run from project root: Rscript 04_Figures/F04/a_script/HRvLR_F04_run.R
#
# Outputs:
#   Main panels A-G : state + within-group volcano-ring contrasts
#   Supp panels A-B : interaction volcano-ring contrasts
#   Main/supp composites via 90_stitch_figure.R

here::i_am(".here")

t0 <- proc.time()

message("\n=== Running panel_A_volcanoes (main + supp) ===")
source(here::here("04_Figures/F04/a_script/panel_A_volcanoes.R"))
message("\n=== Stitching F04 composites ===")
source(here::here("04_Figures/F04/a_script/90_stitch_figure.R"))

elapsed <- proc.time() - t0
message(sprintf("\nAll F04 panels complete (%.0fs)", elapsed["elapsed"]))
