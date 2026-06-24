#!/usr/bin/env Rscript
# F02 global proteome overview — render the individual panels into b_reports.
# No composite; panels render at one shared height and are stitched manually.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

# Wipe stale panel outputs so only the current panel set remains
unlink(setdiff(list.files(RPT_DIR, full.names = TRUE), file.path(RPT_DIR, ".gitkeep")),
  recursive = TRUE
)

source("04_Figures/F02/a_script/panel_A.R") # panel_a_pca
source("04_Figures/F02/a_script/panel_B.R") # panel_b_effect_size
source("04_Figures/F02/a_script/panel_C.R") # panel_c_cv_scatter
source("04_Figures/F02/a_script/panel_F.R") # panel_f_dep_counts
source("04_Figures/F02/a_script/panel_G.R") # panel_g_upset (supp)
source("04_Figures/F02/a_script/panel_H.R") # panel_h_pathways
source("04_Figures/F02/a_script/panel_I.R") # panel_i_venn_training
source("04_Figures/F02/a_script/panel_J.R") # panel_j_venn_acute
source("04_Figures/F02/a_script/panel_M.R") # panel_m_eta2

cat("F02 panels rendered to", RPT_DIR, "\n")
