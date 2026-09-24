# Build phenotype.csv, one row per subject, from the sample sheet. Each trait is computed as
# T2 - T1, except comp_hypertrophy (the source's composite, read from the T2 row) and volume_load.
#
# The MyoVision columns are fibre counts, not areas, although their names say "fCSA": the
# source workbook calls them "Number of fCSA - Mixed (MyoVision)". They run 142-1119 where
# the areas run 3800-10700, and fall as fibre area rises.
#
# ACCUM_VL sits on the T3 row once per subject. It is total kilograms lifted over the
# programme, an exposure rather than an outcome.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
})

meta <- read_csv(here("00_Input", "HRvLR_meta.csv"), show_col_types = FALSE)
trait_cols <- c(
  "fCSA_Type_I_Pre", "fCSA_Type_II_Pre", "fCSA_Mixed_Pre",
  "MyoVision_fCSA_mixed_Pre", "MyoVision_fCSA_Type_I__Pre",
  "mCSA_Pre", "X1RM_Leg_Pre", "X1RM._Ext_Pre"
)
at <- function(timepoint) {
  meta |>
    filter(Timepoint == timepoint) |>
    select(subject = Subject_ID, all_of(trait_cols))
}

deltas <- at("T2") |>
  left_join(at("T1"), by = "subject", suffix = c("_t2", "_t1")) |>
  transmute(
    subject,
    d_fcsa_I = fCSA_Type_I_Pre_t2 - fCSA_Type_I_Pre_t1,
    d_fcsa_II = fCSA_Type_II_Pre_t2 - fCSA_Type_II_Pre_t1,
    d_fcsa_mixed = fCSA_Mixed_Pre_t2 - fCSA_Mixed_Pre_t1,
    d_nfibre_mixed = MyoVision_fCSA_mixed_Pre_t2 - MyoVision_fCSA_mixed_Pre_t1,
    d_nfibre_I = MyoVision_fCSA_Type_I__Pre_t2 - MyoVision_fCSA_Type_I__Pre_t1,
    d_mcsa = mCSA_Pre_t2 - mCSA_Pre_t1,
    d_1rm_legpress = X1RM_Leg_Pre_t2 - X1RM_Leg_Pre_t1,
    d_1rm_ext = X1RM._Ext_Pre_t2 - X1RM._Ext_Pre_t1
  )

phenotype <- meta |>
  distinct(subject = Subject_ID, arm = Group) |>
  left_join(
    meta |>
      filter(Timepoint == "T2") |>
      transmute(subject = Subject_ID, comp_hypertrophy = parse_number(COMP.HYPERTROPHY)),
    by = "subject"
  ) |>
  left_join(deltas, by = "subject") |>
  left_join(
    meta |>
      filter(!is.na(ACCUM_VL)) |>
      transmute(subject = Subject_ID, volume_load = ACCUM_VL),
    by = "subject"
  )

write_csv(phenotype, here("00_Input", "phenotype.csv"))
message("phenotype.csv: ", nrow(phenotype), " subjects x ", ncol(phenotype) - 2, " traits")
