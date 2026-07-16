#!/usr/bin/env Rscript
# Continuous-phenotype prediction from baseline features.
source(here::here("04_Figures", "F06_prediction", "prediction_shared", "a_script", "_leaf.R"))

bundle <- pred_load()
run_cont_leaf(
  bundle, "baseline",
  here::here("04_Figures", "F06_prediction", "prediction_continuous", "baseline")
)
