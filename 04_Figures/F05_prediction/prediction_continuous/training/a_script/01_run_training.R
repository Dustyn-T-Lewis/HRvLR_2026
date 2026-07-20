#!/usr/bin/env Rscript
# Continuous-phenotype prediction from training features.
source(here::here("04_Figures", "F05_prediction", "shared", "a_script", "_leaf.R"))

bundle <- pred_load()
run_cont_leaf(
  bundle, "training",
  here::here("04_Figures", "F05_prediction", "prediction_continuous", "training")
)
