# Run once per proteome clustering (T1, T2, T3, training change, acute change) on that
# clustering's groups C1 and C2: a copy of 02_HRvLR_groups/06_Classification/03_Modules. Each input
# writes to its own folder; stack_inputs() then joins them into c_data/<content>_all_inputs.xlsx
# with a validity tag. How well each module eigengene separates the groups, on the five between-
# group tasks of classify_step().

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ComplexHeatmap)
  library(writexl)
})
source(here("R", "helpers.R"))
source(here("R", "classify.R"))
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

run_input <- function(input) {
  inputs <- c(
    modules = "03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/{input}/modules.rds"
  )
  inputs <- sub("{input}", input, inputs, fixed = TRUE)
  paths <- map_chr(inputs, here)
  stopifnot(file.exists(paths))
  modules <- readRDS(paths[["modules"]])
  values <- modules$eigengenes
  meta <- modules$meta
  stopifnot(identical(colnames(values), meta$sample_id))

  classify_step(
    stage = here("03_Cluster_groups", "Proteome", "06_Classification", "03_Modules", input),
    name = "module_classification", inputs = inputs, paths = paths,
    values = values, meta = meta, feature_label = set_names(rownames(values)),
    feature_noun = "modules", feature_unit = "module", value_label = "eigengene",
    row_groups = list(modules = modules$module_summary$module),
    row_order = "Rows follow module size.",
    groups = c("C1", "C2")
  )
}
purrr::walk(inputs_used, run_input)
stack_inputs(here("03_Cluster_groups", "Proteome", "06_Classification", "03_Modules"))
