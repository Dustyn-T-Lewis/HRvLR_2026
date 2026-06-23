#!/usr/bin/env Rscript
# F02 global proteome overview — render the individual panels (a-e) into b_reports.
# No composite; panels are assembled into the figure manually.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

# Wipe stale panel outputs so only the current panel set remains
unlink(setdiff(list.files(RPT_DIR, full.names = TRUE), file.path(RPT_DIR, ".gitkeep")),
  recursive = TRUE
)

source("04_Figures/F02/a_script/panel_A.R") # panel_a_pca
source("04_Figures/F02/a_script/panel_B.R") # panel_b_effect_size
source("04_Figures/F02/a_script/panel_C.R") # panel_c_cv
source("04_Figures/F02/a_script/panel_D.R") # panel_d_cv_transitions
source("04_Figures/F02/a_script/panel_E.R") # panel_e_intra_variability
source("04_Figures/F02/a_script/panel_F.R") # panel_f_dep_counts
source("04_Figures/F02/a_script/panel_G.R") # panel_g_upset
source("04_Figures/F02/a_script/panel_H.R") # panel_h_pathways
source("04_Figures/F02/a_script/panel_I.R") # panel_i_venn_training
source("04_Figures/F02/a_script/panel_J.R") # panel_j_venn_acute

source("04_Figures/F02/a_script/HRvLR_F02_composite.R") # F02_composite

cat("F02 panels + composite rendered to", RPT_DIR, "\n")
