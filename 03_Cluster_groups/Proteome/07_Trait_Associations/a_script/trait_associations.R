# Each proteome clustering (T1, T2, T3, training, acute), groups C1 and C2: whether a feature's link
# to a trait differs between the groups. The analysis is run_group_trait_associations() in
# 00_Setup/Shared_analyses/run_group_trait_associations.R, shared by every grouping, run once per clustering input;
# stack_inputs() then joins the five workbooks with a validity tag.

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
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

folder <- "03_Cluster_groups/Proteome"
walk(inputs_used, \(input) {
  run_group_trait_associations(
    stage = file.path(folder, "07_Trait_Associations", input),
    design_file = file.path(folder, "02_Design_Contrasts", input, "c_data/design.rds")
  )
})
stack_inputs(here(folder, "07_Trait_Associations"), leave_out = "within_group")
