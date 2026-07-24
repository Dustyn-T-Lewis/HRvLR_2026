suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_rollup.R"))
}))

rollup_root("03_continuous_prediction", "perm_p_q2")
