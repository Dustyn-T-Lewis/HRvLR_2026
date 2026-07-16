#!/usr/bin/env Rscript
# Responder-class (HR vs LR) prediction from baseline features.
source(here::here("04_Figures", "F06_prediction", "prediction_shared", "a_script", "_leaf.R"))

bundle <- pred_load()
run_class_leaf(
  bundle, "baseline",
  here::here("04_Figures", "F06_prediction", "prediction_responder", "baseline")
)
