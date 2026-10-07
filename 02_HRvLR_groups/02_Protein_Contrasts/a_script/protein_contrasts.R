# HR and LR: the nine contrasts on every protein. The analysis is protein_contrasts() in
# R/protein_contrasts.R, shared by every grouping.

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
source(here("R", "protein_contrasts.R"))

protein_contrasts(
  stage = "02_HRvLR_groups/02_Protein_Contrasts",
  design_file = "02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds"
)
