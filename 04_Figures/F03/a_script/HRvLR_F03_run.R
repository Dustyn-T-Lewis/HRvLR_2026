#!/usr/bin/env Rscript
# F03 enrichVolcano ring-volcanoes - compute the fgsea cache (once) and render
# the ring-volcano composites (differential, interaction, responses) into b_reports/.
pacman::p_load(here)
source(here("04_Figures", "F03", "a_script", "HRvLR_F03_setup.R"))
source(here("04_Figures", "F03", "a_script", "panels", "panel_volcanoes.R"))
cat("F03 complete.\n")
