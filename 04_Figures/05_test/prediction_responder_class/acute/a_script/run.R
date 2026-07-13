#!/usr/bin/env Rscript
# Responder-class (HR vs LR) prediction from acute features.
source(here::here("04_Figures", "05_test", "prediction_shared", "_leaf.R"))

bundle <- pred_load()
run_class_leaf(
  bundle, "acute",
  here::here("04_Figures", "05_test", "prediction_responder_class", "acute")
)
