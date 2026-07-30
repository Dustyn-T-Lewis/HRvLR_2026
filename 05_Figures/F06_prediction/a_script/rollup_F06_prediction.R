suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_rollup.R"))
  source(here("functions", "sweep_speccurve.R"))
  source(here("functions", "sweep_manifest.R"))
}))

rollup_root("F06_prediction", "perm_p_q2")
render_root_speccurve("F06_prediction")
build_manifest(
  sweep_root_dir("F06_prediction"), "cont", 792L, "F06_prediction"
)
