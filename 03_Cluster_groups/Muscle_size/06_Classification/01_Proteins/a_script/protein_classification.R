# The muscle-size clusters Hi and Lo: how well each protein separates the groups. The analysis is
# protein_classification() in R/classify.R, shared by every grouping.

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

protein_classification(
  stage = "03_Cluster_groups/Muscle_size/06_Classification/01_Proteins",
  design_file = "03_Cluster_groups/Muscle_size/02_Design_Contrasts/c_data/design.rds"
)
