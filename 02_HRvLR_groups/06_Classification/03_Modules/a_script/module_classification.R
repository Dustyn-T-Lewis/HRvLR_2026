# HR and LR: how well each module eigengene separates the groups. The analysis is
# module_classification() in R/classify.R, shared by every grouping.

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
source(here("R", "helpers.R"))
source(here("R", "classify.R"))

module_classification(
  stage = "02_HRvLR_groups/06_Classification/03_Modules",
  design_file = "02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds"
)
