suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_rollup.R"))
  source(here("functions", "sweep_speccurve.R"))
  source(here("functions", "sweep_manifest.R"))
}))

rollup_root("F05_classification", "perm_p")
render_root_speccurve("F05_classification")
build_manifest(
  sweep_root_dir("F05_classification"), "class", 153L, "F05_classification"
)
