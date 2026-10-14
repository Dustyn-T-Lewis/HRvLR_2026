# The strength clusters Hi and Lo: how well each gene-set score separates the groups. The analysis
# is run_set_discrimination() in 00_Setup/Shared_analyses/run_feature_discrimination.R, shared by every grouping.

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

run_set_discrimination(
  stage = "03_Cluster_groups/Muscle_strength/06_Classification/02_Sets",
  design_file = "03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds"
)
