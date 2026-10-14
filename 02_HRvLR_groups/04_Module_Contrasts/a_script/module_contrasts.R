# HR and LR: the nine contrasts on every module, by eigengene and fry. The analysis is
# run_module_contrasts() in 00_Setup/Shared_analyses/run_module_contrasts.R, shared by every grouping.

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
source(here("00_Setup", "Shared_analyses", "run_module_contrasts.R"))

run_module_contrasts(
  stage = "02_HRvLR_groups/04_Module_Contrasts",
  design_file = "02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds",
  fit_file = "02_HRvLR_groups/02_Protein_Contrasts/c_data/fit.rds"
)
