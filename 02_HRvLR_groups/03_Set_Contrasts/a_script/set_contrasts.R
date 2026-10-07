# HR and LR: the nine contrasts on every gene set, by fgsea and fry. The analysis is set_contrasts()
# in R/set_contrasts.R, shared by every grouping.

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
source(here("R", "helpers.R"))
source(here("R", "set_contrasts.R"))

set_contrasts(
  stage = "02_HRvLR_groups/03_Set_Contrasts",
  design_file = "02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds",
  fit_file = "02_HRvLR_groups/02_Protein_Contrasts/c_data/fit.rds"
)
