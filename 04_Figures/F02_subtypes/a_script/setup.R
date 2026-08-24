# F02 setup: style, the stage-03 clustering tables, and the output paths. Panels
# and 01_run_subtypes.R source this first; idempotent, writes nothing.
#
# Provides: cluster_cells, null_draws, forced_two_group, SPACE_LABELS,
#           F02_RPT, F02_DAT
pacman::p_load(here, dplyr, tibble, readr, openxlsx)

source(here("functions", "shared_style.R"))

F02_RPT <- here("04_Figures", "F02_subtypes", "b_reports")
F02_DAT <- here("04_Figures", "F02_subtypes", "c_data")

SUB_DIR <- here("03_Features", "03_Subtypes", "c_data")

cluster_cells <- as_tibble(read.xlsx(
  file.path(SUB_DIR, "01_subtypes.xlsx"), "cluster_cells"
))
forced_two_group <- as_tibble(read.xlsx(
  file.path(SUB_DIR, "01_subtypes.xlsx"), "forced_two_group"
))
null_draws <- read_csv(file.path(SUB_DIR, "01_null_draws.csv"),
  show_col_types = FALSE
)

SPACE_LABELS <- c(
  eigengenes = "Module eigengenes (12)",
  proteins = "Top-variance proteins (500)"
)

if (!exists("F02_PANELS")) F02_PANELS <- list()
if (!exists("F02_AUDIT")) F02_AUDIT <- list()
