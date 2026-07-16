# F01 setup: style, the phenotype stats, the metadata, and the output paths.
# Panel scripts and 01_run_phenotype.R source this first; it is idempotent. Writes nothing -
# 01_run_phenotype.R creates the output directories.
pacman::p_load(here, dplyr, tidyr, tibble, purrr)

source(here::here("04_Figures", "functions", "shared_style.R"))
source(here::here("04_Figures", "F01_phenotype", "a_script", "phenotype.R"))

F01_RPT <- here::here("04_Figures", "F01_phenotype", "b_reports")
F01_DAT <- here::here("04_Figures", "F01_phenotype", "c_data")

meta <- f01_meta()
if (!exists("F01_PANELS")) F01_PANELS <- list()
if (!exists("F01_AUDIT")) F01_AUDIT <- list()
