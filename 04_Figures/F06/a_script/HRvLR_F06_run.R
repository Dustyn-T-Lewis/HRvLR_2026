#!/usr/bin/env Rscript
# F06 WGCNA module atlas: render panels (A atlas, B ORA, C phenotype) into
# b_reports/panels, then stitch the composite and write the workbook.
pacman::p_load(here)
source(here("04_Figures", "F06", "a_script", "HRvLR_F06_setup.R"))

unlink(
  setdiff(list.files(RPT_DIR, full.names = TRUE), file.path(RPT_DIR, ".gitkeep")),
  recursive = TRUE
)
dir.create(file.path(RPT_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)

source(here("04_Figures", "F06", "a_script", "panels", "panel_A.R"))
source(here("04_Figures", "F06", "a_script", "panels", "panel_B.R"))
source(here("04_Figures", "F06", "a_script", "panels", "panel_C.R"))
source(here("04_Figures", "F06", "a_script", "panels", "panel_S_phenotype_lmm.R"))

source(here("04_Figures", "F06", "a_script", "HRvLR_F06_composite.R"))
cat("F06 complete.\n")
