suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_speccurve.R"))
  source(here("functions", "sweep_manifest.R"))
}))

render_assoc_speccurve("F04_association")
build_manifest(
  sweep_root_dir("F04_association"), "assoc", 420L, "F04_association"
)
