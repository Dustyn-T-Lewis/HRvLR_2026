pacman::p_load(here)
source(here(
  "05_Figures", "F04_association", "a_script", "contrast_level_runner.R"
))

invisible(run_contrast_level("proteins"))
