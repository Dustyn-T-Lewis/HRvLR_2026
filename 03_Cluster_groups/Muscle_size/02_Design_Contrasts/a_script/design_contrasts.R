# The muscle-size clusters Hi and Lo: the design, the two groups and the nine contrasts. The
# analysis is design_contrasts() in R/design_contrasts.R, shared by every grouping.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(purrr)
  library(tibble)
})
source(here("R", "helpers.R"))
source(here("R", "design_contrasts.R"))

design_contrasts(
  stage = "03_Cluster_groups/Muscle_size/02_Design_Contrasts",
  groups = "03_Cluster_groups/Muscle_size/01_Clustering/c_data/groups.rds"
)
