#!/usr/bin/env Rscript
# Responder-class (HR vs LR) prediction from acute features.
source(here::here("04_Figures", "F05_prediction", "shared", "a_script", "_leaf.R"))

bundle <- pred_load()
run_class_leaf(
  bundle, "acute",
  here::here("04_Figures", "F05_prediction", "prediction_responder", "acute")
)
