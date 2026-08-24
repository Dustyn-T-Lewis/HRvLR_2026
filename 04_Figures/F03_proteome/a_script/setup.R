# F03 setup: style, the stage-04 sweep and its permutation check, the stage-05
# pathway calibration, and the output paths. Panels and 01_run_proteome.R source
# this first; idempotent, writes nothing.
#
# Provides: sweep_summary, confirmation, perm_null, fgsea_calibration,
#           fgsea_null, LABEL_NAMES, F03_RPT, F03_DAT
pacman::p_load(here, dplyr, tibble, readr, openxlsx)

source(here("functions", "shared_style.R"))

F03_RPT <- here("04_Figures", "F03_proteome", "b_reports")
F03_DAT <- here("04_Figures", "F03_proteome", "c_data")

PROT_DIR <- here("03_Features", "04_Proteins", "c_data")
PATH_DIR <- here("03_Features", "05_Pathways", "c_data")

sweep_summary <- as_tibble(read.xlsx(
  file.path(PROT_DIR, "01_sweep.xlsx"), "sweep_summary"
))
confirmation <- as_tibble(read.xlsx(
  file.path(PROT_DIR, "02_confirmation.xlsx"), "confirmation"
))
sweep_calibration <- as_tibble(read.xlsx(
  file.path(PROT_DIR, "02_confirmation.xlsx"), "sweep_calibration"
))
perm_null <- read_csv(file.path(PROT_DIR, "02_perm_null_draws.csv"),
  show_col_types = FALSE
)
fgsea_calibration <- as_tibble(read.xlsx(
  file.path(PATH_DIR, "02_fgsea_calibration.xlsx"), "calibration"
))
fgsea_null <- read_csv(file.path(PATH_DIR, "02_fgsea_null_counts.csv"),
  show_col_types = FALSE
)

LABEL_NAMES <- c(
  given = "Given HR/LR", fcsa_I = "fCSA type I", fcsa_II = "fCSA type II",
  fcsa_mixed = "fCSA mixed", nfibre_mixed = "Fibre count, mixed",
  nfibre_I = "Fibre count, type I", mcsa = "Whole-muscle CSA",
  `1rm_legpress` = "1RM leg press", `1rm_ext` = "1RM leg extension",
  volume_load = "Training volume load"
)

# One row per partition, named by the labels that share it, so the panel shows
# 8 tested splits rather than 10 labels of which three are the same split.
partition_names <- function(summary_tbl) {
  summary_tbl |>
    distinct(.data$partition, .data$label) |>
    summarise(
      split = paste(LABEL_NAMES[.data$label], collapse = " / "),
      .by = "partition"
    )
}

if (!exists("F03_PANELS")) F03_PANELS <- list()
if (!exists("F03_AUDIT")) F03_AUDIT <- list()
