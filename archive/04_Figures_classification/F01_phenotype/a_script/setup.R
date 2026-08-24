# F01 setup: style, the stage-01 audit tables, and the output paths. Panels and
# 01_run_phenotype.R source this first; it is idempotent and writes nothing.
#
# Provides: pheno, change_summary, label_separation, composite_structure,
#           composite_modality, ranked, TRAIT_LABELS, F01_RPT, F01_DAT
pacman::p_load(here, dplyr, tibble, readr, openxlsx)

source(here("functions", "shared_style.R"))

F01_RPT <- here("04_Figures", "F01_phenotype", "b_reports")
F01_DAT <- here("04_Figures", "F01_phenotype", "c_data")

AUDIT_BOOK <- here(
  "03_Features", "01_Responsiveness", "c_data", "01_label_audit.xlsx"
)

pheno <- read_csv(here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
)

change_summary <- as_tibble(read.xlsx(AUDIT_BOOK, "change_summary"))
label_separation <- as_tibble(read.xlsx(AUDIT_BOOK, "label_separation"))
composite_structure <- as_tibble(read.xlsx(AUDIT_BOOK, "composite_structure"))
composite_modality <- as_tibble(read.xlsx(AUDIT_BOOK, "composite_modality"))

ranked <- pheno |>
  arrange(.data$comp_hypertrophy) |>
  mutate(subject = forcats::fct_inorder(.data$subject))

TRAIT_LABELS <- c(
  comp_hypertrophy = "Composite hypertrophy",
  d_fcsa_I = "fCSA type I",
  d_fcsa_II = "fCSA type II",
  d_fcsa_mixed = "fCSA mixed",
  d_nfibre_mixed = "Fibre count, mixed",
  d_nfibre_I = "Fibre count, type I",
  d_mcsa = "Whole-muscle CSA",
  d_1rm_legpress = "1RM leg press",
  d_1rm_ext = "1RM leg extension",
  volume_load = "Training volume load"
)

# Both effect-size panels use this one order so a trait keeps its row between
# them and the reversal between "changed" and "separated" reads as a flip in
# place rather than a reshuffle.
# volume_load is a total, not a change, so it has no row in change_summary and
# is appended at the end rather than ordered among the change scores.
TRAIT_ORDER <- change_summary |>
  filter(.data$trait != "comp_hypertrophy") |>
  arrange(.data$mean_d) |>
  pull(.data$trait) |>
  (\(x) unname(TRAIT_LABELS[c(x, "volume_load")]))()

if (!exists("F01_PANELS")) F01_PANELS <- list()
if (!exists("F01_AUDIT")) F01_AUDIT <- list()
