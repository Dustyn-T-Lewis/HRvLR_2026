# Each proteome clustering (T1, T2, T3, training, acute), groups C1 and C2: the design, the two
# groups and the nine contrasts. The analysis is design_contrasts() in R/design_contrasts.R, shared
# by every grouping, run once per clustering input; stack_inputs() then joins the five workbooks
# with a validity tag.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(purrr)
  library(tibble)
})
source(here("R", "helpers.R"))
source(here("R", "design_contrasts.R"))
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

folder <- "03_Cluster_groups/Proteome"
walk(inputs_used, \(input) {
  design_contrasts(
    stage = file.path(folder, "02_Design_Contrasts", input),
    groups = file.path(folder, "01_Clustering/c_data/groups", paste0(input, ".rds"))
  )
})
stack_inputs(here(folder, "02_Design_Contrasts"))
