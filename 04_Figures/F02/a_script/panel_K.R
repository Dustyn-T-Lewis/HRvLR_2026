# F02 Panel K (pilot): HR vs LR training-response overlap at nominal p < 0.05 ----
# Looser-threshold companion to panel_I (Pi < 0.05); supplementary.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")
if (!exists("venn_panel")) source("04_Figures/F02/a_script/panel_I.R")

venn_panel("Training_HR", "Training_LR",
  "Training Response (T2 - T1, p < 0.05)",
  metric = "p",
  file.path(RPT_DIR, "panel_k_venn_training_p")
)
cat("F02 Panel K done.\n")
