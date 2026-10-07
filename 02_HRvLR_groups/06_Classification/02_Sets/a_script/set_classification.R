# HR and LR: how well each gene-set score separates the groups. The analysis is set_classification()
# in R/classify.R, shared by every grouping.

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

set_classification(
  stage = "02_HRvLR_groups/06_Classification/02_Sets",
  design_file = "02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds"
)
