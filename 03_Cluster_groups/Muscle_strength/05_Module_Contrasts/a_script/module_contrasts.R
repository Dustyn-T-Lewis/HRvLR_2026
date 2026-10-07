# The strength clusters Hi and Lo: the nine contrasts on every module, by eigengene and fry. The
# analysis is module_contrasts() in R/module_contrasts.R, shared by every grouping.

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
source(here("R", "module_contrasts.R"))

module_contrasts(
  stage = "03_Cluster_groups/Muscle_strength/05_Module_Contrasts",
  design_file = "03_Cluster_groups/Muscle_strength/02_Design_Contrasts/c_data/design.rds",
  fit_file = "03_Cluster_groups/Muscle_strength/03_Protein_Contrasts/c_data/fit.rds"
)
