# How well each module eigengene separates the groups, on the five between-group tasks of
# classify_step().

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

inputs <- c(modules = "01_Trunk/06_Module_Building/c_data/modules.rds")
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
modules <- readRDS(paths[["modules"]])
values <- modules$eigengenes
meta <- modules$meta
stopifnot(identical(colnames(values), meta$sample_id))

classify_step(
  stage = here("02_HRvLR_groups", "06_Classification", "03_Modules"),
  name = "module_classification", inputs = inputs, paths = paths,
  values = values, meta = meta, feature_label = set_names(rownames(values)),
  feature_noun = "modules", feature_unit = "module", value_label = "eigengene",
  row_groups = list(modules = modules$module_summary$module), row_order = "Rows follow module size."
)
