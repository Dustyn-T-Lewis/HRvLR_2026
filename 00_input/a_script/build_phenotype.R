pacman::p_load(dplyr, readr)

build_phenotype_table <- function(meta_path) {
  meta <- read_csv(meta_path, show_col_types = FALSE)

  trait_cols <- c(
    "fCSA_Type_I_Pre", "fCSA_Type_II_Pre",
    "mCSA_Pre", "X1RM_Leg_Pre", "X1RM._Ext_Pre"
  )

  arm <- meta |>
    distinct(Subject_ID, Group) |>
    rename(subject = Subject_ID, group_arm = Group)

  comp <- meta |>
    filter(Timepoint == "T2") |>
    transmute(subject = Subject_ID, comp_hypertrophy = parse_number(COMP.HYPERTROPHY))

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
      d_mcsa = mCSA_Pre_t2 - mCSA_Pre_t1,
      d_1rm_legpress = X1RM_Leg_Pre_t2 - X1RM_Leg_Pre_t1,
      d_1rm_ext = `X1RM._Ext_Pre_t2` - `X1RM._Ext_Pre_t1`
    ) |>
    select(subject, d_fcsa_I, d_fcsa_II, d_mcsa, d_1rm_legpress, d_1rm_ext)

  arm |>
    left_join(comp, by = "subject") |>
    left_join(deltas, by = "subject")
}
