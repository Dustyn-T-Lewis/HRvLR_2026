#!/usr/bin/env Rscript
# Responder-class (HR vs LR) prediction from training features.
source(here::here("04_Figures", "F06_prediction", "prediction_shared", "_leaf.R"))

bundle <- pred_load()
run_class_leaf(
  bundle, "training",
  here::here("04_Figures", "F06_prediction", "prediction_responder", "training")
)
