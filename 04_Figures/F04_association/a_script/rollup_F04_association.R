suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_manifest.R"))
}))

build_manifest(
  sweep_root_dir("F04_association"), "assoc", 420L, "F04_association"
)
