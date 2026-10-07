# The strength clusters Hi and Lo: how well each module eigengene separates the groups. The analysis
# is module_classification() in R/classify.R, shared by every grouping.

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
  stage = "03_Cluster_groups/Muscle_strength/06_Classification/03_Modules",
  design_file = "03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds"
)
