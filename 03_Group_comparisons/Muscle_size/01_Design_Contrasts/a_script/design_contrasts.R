# The muscle-size clusters Hi and Lo: the design, the two groups and the nine contrasts. The
# analysis is run_design_contrasts() in 00_Setup/Shared_analyses/run_design_contrasts.R, shared by every grouping.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(purrr)
  library(tibble)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_analyses", "run_design_contrasts.R"))

run_design_contrasts(
  stage = "03_Group_comparisons/Muscle_size/01_Design_Contrasts",
  groups = "02_Groups/Muscle_size/c_data/groups.rds"
)
