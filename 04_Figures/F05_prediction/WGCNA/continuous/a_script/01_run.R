# Continuous-phenotype prediction from the WGCNA module eigengenes, across phases.
source(here::here("04_Figures", "F05_prediction", "shared", "a_script", "_leaf.R"))
run_cont_arm(
  pred_load(),
  here::here("04_Figures", "F05_prediction", "WGCNA", "continuous"),
  space = "eigengenes"
)
