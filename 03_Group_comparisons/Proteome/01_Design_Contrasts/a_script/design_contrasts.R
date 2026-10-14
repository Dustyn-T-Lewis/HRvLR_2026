# Each proteome clustering (T1, T2, T3, training, acute), groups C1 and C2: the design, the two
# groups and the nine contrasts. The analysis is run_design_contrasts() in 00_Setup/Shared_analyses/run_design_contrasts.R, shared
# by every grouping, run once per clustering input; stack_inputs() then joins the five workbooks
# with a validity tag.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(purrr)
  library(tibble)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_analyses", "run_design_contrasts.R"))
source(here("00_Setup", "Shared_functions", "proteome_inputs.R"))

folder <- "03_Group_comparisons/Proteome"
walk(inputs_used, \(input) {
  run_design_contrasts(
    stage = file.path(folder, "01_Design_Contrasts", input),
    groups = file.path("02_Groups/Proteome/c_data/groups", paste0(input, ".rds"))
  )
})
stack_inputs(here(folder, "01_Design_Contrasts"))
