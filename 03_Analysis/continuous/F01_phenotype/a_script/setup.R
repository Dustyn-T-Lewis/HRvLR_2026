# F01 setup, continuous tree: style, the phenotype stats, the metadata, and
# the output paths. Reuses the same phenotype_helpers.R as the categorical
# tree -- the composite score and the raw outcome table are the same data
# regardless of tree; only the group-contrast functions in that file go
# unused here.
pacman::p_load(here, dplyr, tidyr, tibble, purrr)

source(here::here("functions", "shared_style.R"))
source(here::here("functions", "phenotype_helpers.R"))

F01_RPT <- here::here("03_Analysis", "continuous", "F01_phenotype", "b_reports")
F01_DAT <- here::here("03_Analysis", "continuous", "F01_phenotype", "c_data")

meta <- f01_meta()
if (!exists("F01_PANELS")) F01_PANELS <- list()
if (!exists("F01_AUDIT")) F01_AUDIT <- list()
