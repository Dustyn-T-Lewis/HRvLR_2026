# The strength clusters Hi and Lo: how well each protein separates the groups. The analysis is
# run_protein_discrimination() in 00_Setup/Shared_analyses/run_feature_discrimination.R, shared by every grouping.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ComplexHeatmap)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_analyses", "run_feature_discrimination.R"))

run_protein_discrimination(
  stage = "03_Group_comparisons/Muscle_strength/05_Classification/01_Proteins",
  design_file = "03_Group_comparisons/Muscle_strength/01_Design_Contrasts/c_data/design.rds"
)
