#!/usr/bin/env Rscript
# HRvLR_F04_run.R — render the nine volcano-ring panels (a-i) into b_reports.
# No composite; panels are assembled into the figure manually.
# Run from project root: Rscript 04_Figures/F04/a_script/HRvLR_F04_run.R

here::i_am(".here")

t0 <- proc.time()
source(here::here("04_Figures/F04/a_script/panel_A_volcanoes.R"))
elapsed <- proc.time() - t0
message(sprintf("All F04 panels complete (%.0fs)", elapsed["elapsed"]))
