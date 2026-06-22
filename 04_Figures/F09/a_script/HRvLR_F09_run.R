#!/usr/bin/env Rscript
# F09 WGCNA module figure — render the standalone panels into b_reports.
# No composite; panels are assembled into the figure manually.
# Requires the WGCNA artifacts from 04_Figures/F06/a_script/HRvLR_WGCNA_run.R.

setwd(rprojroot::find_rstudio_root_file())

RPT <- "04_Figures/F09/b_reports"
unlink(setdiff(list.files(RPT, full.names = TRUE), file.path(RPT, ".gitkeep")),
  recursive = TRUE
)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

source("04_Figures/F09/a_script/panel_A.R") # panel_a_module_trait
source("04_Figures/F09/a_script/panel_B_triptych.R") # panel_b_triptych
source("04_Figures/F09/a_script/panel_D.R") # panel_c_hub_<module>

cat("F09 WGCNA panels rendered to", RPT, "\n")
