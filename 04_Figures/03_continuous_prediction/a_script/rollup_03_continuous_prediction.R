suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_rollup.R"))
  source(here("functions", "sweep_speccurve.R"))
}))

rollup_root("03_continuous_prediction", "perm_p_q2")
render_root_speccurve("03_continuous_prediction")
