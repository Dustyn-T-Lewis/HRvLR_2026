# F01 setup: load builders, the metadata, and the output paths. Panel scripts and
# run.R source this first; it is idempotent.
source(here::here("04_Figures", "F01", "a_script", "functions", "build.R"))

F01_RPT <- here::here("04_Figures", "F01", "b_reports")
F01_DAT <- here::here("04_Figures", "F01", "c_data")
dir.create(file.path(F01_RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(F01_RPT, "supp"), recursive = TRUE, showWarnings = FALSE)
dir.create(F01_DAT, recursive = TRUE, showWarnings = FALSE)

meta <- f01_meta()
if (!exists("F01_PANELS")) F01_PANELS <- list()
if (!exists("F01_AUDIT")) F01_AUDIT <- list()
