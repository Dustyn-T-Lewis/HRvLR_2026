# Each proteome clustering (T1, T2, T3, training, acute), groups C1 and C2: how well each module
# eigengene separates the groups. The analysis is run_module_discrimination() in 00_Setup/Shared_analyses/run_feature_discrimination.R, shared by
# every grouping, run once per clustering input; stack_inputs() then joins the five workbooks with a
# validity tag.

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
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_analyses", "run_feature_discrimination.R"))
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

folder <- "03_Cluster_groups/Proteome"
walk(inputs_used, \(input) {
  run_module_discrimination(
    stage = file.path(folder, "06_Classification/03_Modules", input),
    design_file = file.path(folder, "02_Design_Contrasts", input, "c_data/design.rds")
  )
})
stack_inputs(here(folder, "06_Classification/03_Modules"))
