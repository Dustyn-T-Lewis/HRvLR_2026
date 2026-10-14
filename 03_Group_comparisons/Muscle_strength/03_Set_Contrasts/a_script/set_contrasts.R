# The strength clusters Hi and Lo: the nine contrasts on every gene set, by fgsea and fry. The
# analysis is run_set_contrasts() in 00_Setup/Shared_analyses/run_set_contrasts.R, shared by every grouping.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ComplexHeatmap)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_analyses", "run_set_contrasts.R"))

run_set_contrasts(
  stage = "03_Group_comparisons/Muscle_strength/03_Set_Contrasts",
  design_file = "03_Group_comparisons/Muscle_strength/01_Design_Contrasts/c_data/design.rds",
  fit_file = "03_Group_comparisons/Muscle_strength/02_Protein_Contrasts/c_data/fit.rds"
)
