#!/usr/bin/env Rscript
# Responder-class (HR vs LR) prediction from training features.
source(here::here("05_test", "prediction_shared", "_leaf.R"))

bundle <- pred_load()
run_class_leaf(
  bundle, "training",
  here::here("05_test", "prediction_responder_class", "training")
)
