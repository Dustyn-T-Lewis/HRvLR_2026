# F03 setup: style, the module atlas, the association sweep restricted to
# modules, and its calibration. Panels source this first; writes nothing.
#
# Provides: atlas, module_results, confirmation, sweep_calibration, robustness,
#           TRAIT_LABELS, WINDOW_LABELS, F03_RPT, F03_DAT
pacman::p_load(here, dplyr, tibble, readr, openxlsx)

source(here("functions", "shared_style.R"))
source(here("functions", "association.R"))

F03_RPT <- here("04_Figures", "F03_modules", "b_reports")
F03_DAT <- here("04_Figures", "F03_modules", "c_data")

FEAT <- here("03_Features", "c_data")

atlas <- as_tibble(read.xlsx(
  file.path(FEAT, "05_module_annotation.xlsx"), "module_atlas"
))
module_results <- read_csv(file.path(FEAT, "02_association_full.csv"),
  show_col_types = FALSE
) |>
  filter(.data$level == "modules")
confirmation <- read_csv(file.path(FEAT, "03_confirmation.csv"),
  show_col_types = FALSE
) |>
  filter(.data$level == "modules")
sweep_calibration <- as_tibble(read.xlsx(
  file.path(FEAT, "03_confirmation.xlsx"), "sweep_calibration"
))
robustness <- as_tibble(read.xlsx(
  file.path(FEAT, "04_hit_robustness.xlsx"), "robustness"
))

TRAIT_LABELS <- c(
  comp_hypertrophy = "Composite hypertrophy",
  d_fcsa_I = "fCSA type I", d_fcsa_II = "fCSA type II",
  d_fcsa_mixed = "fCSA mixed", d_nfibre_mixed = "Fibre count, mixed",
  d_nfibre_I = "Fibre count, type I", d_mcsa = "Whole-muscle CSA",
  d_1rm_legpress = "1RM leg press", d_1rm_ext = "1RM leg extension",
  volume_load = "Training volume load"
)

WINDOW_LABELS <- c(
  T1 = "T1", T2 = "T2", T3 = "T3",
  training = "Training", acute = "Acute", total = "Total"
)

# Short names from the GO cellular-component term, which is the sub-ontology
# that named these modules most sharply.
module_label <- function(mod) {
  cc <- atlas$CC[match(mod, atlas$module)]
  bp <- atlas$BP[match(mod, atlas$module)]
  term <- ifelse(is.na(cc), bp, cc)
  term <- sub(" \\(q=.*$", "", term)
  ifelse(is.na(term), mod, paste0(mod, ": ", term))
}

if (!exists("F03_PANELS")) F03_PANELS <- list()
if (!exists("F03_AUDIT")) F03_AUDIT <- list()
