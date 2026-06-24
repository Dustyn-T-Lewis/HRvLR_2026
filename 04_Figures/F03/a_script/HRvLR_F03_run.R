#!/usr/bin/env Rscript
# F03 module-phenotype linkage — render the individual panels into b_reports.
# No composite; panels render at one shared height and are stitched manually.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F03/a_script/HRvLR_F03_setup.R")

# Wipe stale panel outputs so only the current panel set remains
unlink(
  setdiff(
    list.files("04_Figures/F03/b_reports", full.names = TRUE),
    file.path("04_Figures/F03/b_reports", ".gitkeep")
  ),
  recursive = TRUE
)

source("04_Figures/F03/a_script/panel_A.R") # panel_a_wgcna_map
source("04_Figures/F03/a_script/panel_B.R") # panel_b_trajectories
source("04_Figures/F03/a_script/panel_C.R") # panel_c_mfuzz_pca
source("04_Figures/F03/a_script/panel_D.R") # panel_d_ora
source("04_Figures/F03/a_script/panel_E.R") # panel_e_verdict

cat("F03 complete.\n")
