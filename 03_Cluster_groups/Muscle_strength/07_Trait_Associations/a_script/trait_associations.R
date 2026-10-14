# The strength clusters Hi and Lo: whether a feature's link to a trait differs between the groups.
# The analysis is run_group_trait_associations() in 00_Setup/Shared_analyses/run_group_trait_associations.R, shared by every
# grouping.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_analyses", "run_group_trait_associations.R"))

run_group_trait_associations(
  stage = "03_Cluster_groups/Muscle_strength/07_Trait_Associations",
  design_file = "03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds"
)
