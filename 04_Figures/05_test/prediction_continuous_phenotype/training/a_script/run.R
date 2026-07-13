#!/usr/bin/env Rscript
# Continuous-phenotype prediction from training features.
source(here::here("04_Figures", "05_test", "prediction_shared", "_leaf.R"))

bundle <- pred_load()
run_cont_leaf(
  bundle, "training",
  here::here("04_Figures", "05_test", "prediction_continuous_phenotype", "training")
)
