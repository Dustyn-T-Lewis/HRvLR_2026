# Responder-class prediction from singscore pathway scores, across phases.
source(here::here("04_Figures", "F05_prediction", "shared", "a_script", "_leaf.R"))
run_class_arm(
  pred_load(),
  here::here("04_Figures", "F05_prediction", "singscore", "responder"),
  space = "singscore"
)
