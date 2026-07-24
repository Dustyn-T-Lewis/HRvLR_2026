suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_rollup.R"))
}))

rollup_root("02_classification", "perm_p")
