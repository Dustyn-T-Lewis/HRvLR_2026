# Baseline continuous prediction of adaptation deltas (spls engine).
source(here::here("functions", "pred_leaf.R"))
run_cont_leaf(
  pred_load(),
  model = "spls",
  out_dir = here::here(
    "04_Figures", "F05_prediction", "phenotype", "baseline", "spls"
  )
)
