# HR and LR: the nine contrasts on every protein. The analysis is run_protein_contrasts() in
# 00_Setup/Shared_analyses/run_protein_contrasts.R, shared by every grouping.

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
source(here("00_Setup", "Shared_analyses", "run_protein_contrasts.R"))

run_protein_contrasts(
  stage = "03_Group_comparisons/HRvLR/02_Protein_Contrasts",
  design_file = "03_Group_comparisons/HRvLR/01_Design_Contrasts/c_data/design.rds"
)
