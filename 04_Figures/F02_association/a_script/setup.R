# F02 setup: style, the association sweep and its calibration, the output
# paths. Panels and 01_run_association.R source this first; writes nothing.
#
# Provides: summary_tbl, survivors, results, confirmation, sweep_calibration,
#           pheno, TRAIT_LABELS, LEVEL_LABELS, F02_RPT, F02_DAT
pacman::p_load(here, dplyr, tibble, readr, openxlsx)

source(here("functions", "shared_style.R"))
source(here("functions", "association.R"))

F02_RPT <- here("04_Figures", "F02_association", "b_reports")
F02_DAT <- here("04_Figures", "F02_association", "c_data")

FEAT_DIR <- here("03_Features", "c_data")

summary_tbl <- as_tibble(read.xlsx(
  file.path(FEAT_DIR, "02_association.xlsx"), "association_summary"
))
survivors <- as_tibble(read.xlsx(
  file.path(FEAT_DIR, "02_association.xlsx"), "survivors"
))
results <- read_csv(file.path(FEAT_DIR, "02_association_full.csv"),
  show_col_types = FALSE
)

CONFIRM_BOOK <- file.path(FEAT_DIR, "03_confirmation.xlsx")
confirmation <- if (file.exists(CONFIRM_BOOK)) {
  as_tibble(read.xlsx(CONFIRM_BOOK, "confirmation"))
} else {
  NULL
}
sweep_calibration <- if (file.exists(CONFIRM_BOOK)) {
  as_tibble(read.xlsx(CONFIRM_BOOK, "sweep_calibration"))
} else {
  NULL
}

pheno <- phenotype_table()

TRAIT_LABELS <- c(
  comp_hypertrophy = "Composite hypertrophy",
  d_fcsa_I = "fCSA type I", d_fcsa_II = "fCSA type II",
  d_fcsa_mixed = "fCSA mixed", d_nfibre_mixed = "Fibre count, mixed",
  d_nfibre_I = "Fibre count, type I", d_mcsa = "Whole-muscle CSA",
  d_1rm_legpress = "1RM leg press", d_1rm_ext = "1RM leg extension",
  volume_load = "Training volume load"
)

LEVEL_LABELS <- c(
  proteins = "Proteins (1900)", modules = "Modules (12)",
  pathways = "Pathways (57)"
)

WINDOW_LABELS <- c(
  T1 = "Level at T1", T2 = "Level at T2", T3 = "Level at T3",
  training = "Training (T2-T1)", acute = "Acute (T3-T2)",
  total = "Total (T3-T1)"
)

# Level windows carry a value, change windows a difference; panels that name
# the y axis need to say which.
window_family <- function(w) ifelse(w %in% LEVEL_WINDOWS, "level", "change")

if (!exists("F02_PANELS")) F02_PANELS <- list()
if (!exists("F02_AUDIT")) F02_AUDIT <- list()
