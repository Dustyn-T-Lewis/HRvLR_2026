#!/usr/bin/env Rscript
# F05 Acute concordance/divergence - render the panels, then stitch + workbook.
pacman::p_load(here)
source(here("04_Figures", "F05", "a_script", "HRvLR_F05_setup.R"))

source(here("04_Figures", "F05", "a_script", "panels", "panel_A_quadrant_ora.R"))
source(here("04_Figures", "F05", "a_script", "panels", "panel_B_nes_scatter.R"))

source(here("04_Figures", "F05", "a_script", "supp", "supp_fry.R"))
source(here("04_Figures", "F05", "a_script", "supp", "supp_rrho2.R"))
source(here("04_Figures", "F05", "a_script", "supp", "supp_threshold.R"))
source(here("04_Figures", "F05", "a_script", "supp", "supp_bootstrap.R"))

source(here("04_Figures", "F05", "a_script", "HRvLR_F05_composite.R"))
cat("F05 complete.\n")
