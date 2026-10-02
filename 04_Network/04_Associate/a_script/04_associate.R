# Ties each module eigengene to phenotype by the sample model and the change scores that
# R/associate.R defines.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(writexl)
  library(grid)
  library(ComplexHeatmap)
})
source(here("R", "helpers.R"))
source(here("R", "associate.R"))

inputs <- c(
  modules = "04_Network/01_Build_Modules/c_data/modules.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
modules <- readRDS(paths[["modules"]])
values <- modules$eigengenes
meta <- as_tibble(modules$meta)
stopifnot(identical(colnames(values), meta$sample_id))

associate_step(
  stage = here("04_Network", "04_Associate"), inputs = inputs, paths = paths,
  values = values, meta = meta, feature_label = set_names(rownames(values)),
  feature_noun = "modules", feature_unit = "module", value_label = "eigengene",
  row_groups = list(modules = modules$module_summary$module),
  row_order = "Rows follow module size.",
  all_arm_rows = TRUE
)
