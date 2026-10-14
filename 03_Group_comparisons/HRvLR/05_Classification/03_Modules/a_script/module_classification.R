# HR and LR: how well each module eigengene separates the groups. The analysis is
# run_module_discrimination() in 00_Setup/Shared_analyses/run_feature_discrimination.R, shared by every grouping.

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

run_module_discrimination(
  stage = "03_Group_comparisons/HRvLR/05_Classification/03_Modules",
  design_file = "03_Group_comparisons/HRvLR/01_Design_Contrasts/c_data/design.rds"
)
