# F02 Panel I: HR vs LR training-response overlap, baseline-anchored (3-set Euler) ----
# Shared vs responder-specific proteins changing with training (T2 - T1, Pi < 0.05),
# overlaid on the baseline HR-vs-LR difference set, with a companion up/down bar.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")
if (!exists("venn3_panel")) source("04_Figures/F02/a_script/functions/venn.R")

venn3_panel(
  "Training_HR", "Training_LR", "Baseline_HRvLR",
  "Training Response (T2 - T1)",
  file.path(RPT_DIR, "panel_i_venn_training")
)
cat("F02 Panel I done.\n")
