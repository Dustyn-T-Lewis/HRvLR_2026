#!/usr/bin/env Rscript
# Responder-class (HR vs LR) prediction from baseline features.
source(here::here("05_test", "prediction_shared", "_leaf.R"))

bundle <- pred_load()
run_class_leaf(
  bundle, "baseline",
  here::here("05_test", "prediction_responder_class", "baseline")
)
