pacman::p_load(dplyr, readr)

# The MyoVision columns are fibre counts, not areas, despite carrying "fCSA" in
# their meta names: the source workbook calls them "Number of fCSA - Mixed
# (MyoVision)". That is why they run 142-1119 where the areas run 3800-10700,
# and why they correlate negatively with area - larger fibres, fewer of them in
# the imaged field.
#
# ACCUM_VL is recorded once per subject on the T3 row and is an exposure, not an
# outcome: total kg lifted over the programme, spanning 3.5-fold across the
# cohort. It is the only variable here that describes what a subject did rather
# than what happened to them.

build_phenotype_table <- function(meta_path) {
  meta <- read_csv(meta_path, show_col_types = FALSE)

  trait_cols <- c(
    "fCSA_Type_I_Pre", "fCSA_Type_II_Pre", "fCSA_Mixed_Pre",
    "MyoVision_fCSA_mixed_Pre", "MyoVision_fCSA_Type_I__Pre",
    "mCSA_Pre", "X1RM_Leg_Pre", "X1RM._Ext_Pre"
  )

  arm <- meta |>
    distinct(Subject_ID, Group) |>
    rename(subject = Subject_ID, group_arm = Group)

  comp <- meta |>
    filter(Timepoint == "T2") |>
    transmute(
      subject = Subject_ID,
      comp_hypertrophy = parse_number(COMP.HYPERTROPHY)
    )

  volume <- meta |>
    filter(!is.na(ACCUM_VL)) |>
    transmute(subject = Subject_ID, volume_load = ACCUM_VL)

  t1 <- meta |>
    filter(Timepoint == "T1") |>
    select(subject = Subject_ID, all_of(trait_cols))
  t2 <- meta |>
    filter(Timepoint == "T2") |>
    select(subject = Subject_ID, all_of(trait_cols))

  deltas <- t2 |>
    left_join(t1, by = "subject", suffix = c("_t2", "_t1")) |>
    mutate(
      d_fcsa_I = fCSA_Type_I_Pre_t2 - fCSA_Type_I_Pre_t1,
      d_fcsa_II = fCSA_Type_II_Pre_t2 - fCSA_Type_II_Pre_t1,
      d_fcsa_mixed = fCSA_Mixed_Pre_t2 - fCSA_Mixed_Pre_t1,
      d_nfibre_mixed =
        MyoVision_fCSA_mixed_Pre_t2 - MyoVision_fCSA_mixed_Pre_t1,
      d_nfibre_I =
        `MyoVision_fCSA_Type_I__Pre_t2` - `MyoVision_fCSA_Type_I__Pre_t1`,
      d_mcsa = mCSA_Pre_t2 - mCSA_Pre_t1,
      d_1rm_legpress = X1RM_Leg_Pre_t2 - X1RM_Leg_Pre_t1,
      d_1rm_ext = `X1RM._Ext_Pre_t2` - `X1RM._Ext_Pre_t1`
    ) |>
    select(
      subject, d_fcsa_I, d_fcsa_II, d_fcsa_mixed, d_nfibre_mixed, d_nfibre_I,
      d_mcsa, d_1rm_legpress, d_1rm_ext
    )

  arm |>
    left_join(comp, by = "subject") |>
    left_join(deltas, by = "subject") |>
    left_join(volume, by = "subject")
}
