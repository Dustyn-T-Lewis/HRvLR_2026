suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_rollup.R"))
  source(here("functions", "sweep_speccurve.R"))
}))

rollup_root("02_classification", "perm_p")
render_root_speccurve("02_classification")
