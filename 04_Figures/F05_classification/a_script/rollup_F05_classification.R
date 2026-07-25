suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_rollup.R"))
  source(here("functions", "sweep_speccurve.R"))
}))

rollup_root("F05_classification", "perm_p")
render_root_speccurve("F05_classification")
