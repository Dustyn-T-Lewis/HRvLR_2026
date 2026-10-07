# Each proteome clustering (T1, T2, T3, training, acute), groups C1 and C2: the nine contrasts on
# every gene set, by fgsea and fry. The analysis is set_contrasts() in R/set_contrasts.R, shared by
# every grouping, run once per clustering input; stack_inputs() then joins the five workbooks with a
# validity tag.

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
source(here("R", "set_contrasts.R"))
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

folder <- "03_Cluster_groups/Proteome"
walk(inputs_used, \(input) {
  set_contrasts(
    stage = file.path(folder, "04_Set_Contrasts", input),
    design_file = file.path(folder, "02_Design_Contrasts", input, "c_data/design.rds"),
    fit_file = file.path(folder, "03_Protein_Contrasts", input, "c_data/fit.rds")
  )
})
stack_inputs(here(folder, "04_Set_Contrasts"))
