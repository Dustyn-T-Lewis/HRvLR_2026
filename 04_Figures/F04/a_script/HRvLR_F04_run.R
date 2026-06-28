#!/usr/bin/env Rscript
# F04 Training concordance/divergence - render the panels, then stitch + workbook.
pacman::p_load(here)
source(here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))

source(here("04_Figures", "F04", "a_script", "panels", "panel_A_quadrant_ora.R"))
source(here("04_Figures", "F04", "a_script", "panels", "panel_B_nes_scatter.R"))

source(here("04_Figures", "F04", "a_script", "supp", "supp_fry.R"))
source(here("04_Figures", "F04", "a_script", "supp", "supp_rrho2.R"))
source(here("04_Figures", "F04", "a_script", "supp", "supp_threshold.R"))
source(here("04_Figures", "F04", "a_script", "supp", "supp_bootstrap.R"))

source(here("04_Figures", "F04", "a_script", "HRvLR_F04_composite.R"))
cat("F04 complete.\n")
