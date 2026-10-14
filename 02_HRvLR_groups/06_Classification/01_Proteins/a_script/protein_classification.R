# HR and LR: how well each protein separates the groups. The analysis is run_protein_discrimination() in
# 00_Setup/Shared_analyses/run_feature_discrimination.R, shared by every grouping.

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
  stage = "02_HRvLR_groups/06_Classification/01_Proteins",
  design_file = "02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds"
)
