# The strength clusters Hi and Lo: how well each gene-set score separates the groups. The analysis
# is set_classification() in R/classify.R, shared by every grouping.

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
  stage = "03_Cluster_groups/Muscle_strength/06_Classification/02_Sets",
  design_file = "03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds"
)
