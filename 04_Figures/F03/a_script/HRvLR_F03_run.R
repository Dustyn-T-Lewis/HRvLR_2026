#!/usr/bin/env Rscript
# F03 enrichVolcano ring-volcanoes - compute the fgsea cache (once) and render
# the 9-contrast ring-volcano grid into b_reports/panels.
pacman::p_load(here)
source(here("04_Figures", "F03", "a_script", "HRvLR_F03_setup.R"))
source(here("04_Figures", "F03", "a_script", "panels", "panel_volcanoes.R"))
cat("F03 complete.\n")
