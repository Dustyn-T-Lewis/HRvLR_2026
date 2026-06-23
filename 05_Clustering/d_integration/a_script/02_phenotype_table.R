# Per-subject phenotype table for module–trait linkage.
# Outputs: 05_Clustering/d_integration/c_data/phenotype.csv

pacman::p_load(here, dplyr, tidyr, readr)

meta <- read_csv(
  here("00_input", "HRvLR_meta.csv"),
  show_col_types = FALSE
)

trait_cols <- c(
  "fCSA_Type_I_Pre",
  "fCSA_Type_II_Pre",
  "MyoVision_fCSA_Type_I__Pre",
  "mCSA_Pre",
  "X1RM_Leg_Pre",
  "X1RM._Ext_Pre"
)

delta_names <- c(
  "d_fcsa_I",
  "d_fcsa_II",
  "d_myovision_fcsa_I",
  "d_mcsa",
  "d_1rm_legpress",
  "d_1rm_ext"
)

arm <- meta |>
  distinct(Subject_ID, Group) |>
  rename(subject = Subject_ID, group_arm = Group)

comp <- meta |>
  filter(Timepoint == "T2") |>
  select(Subject_ID, COMP.HYPERTROPHY) |>
  mutate(comp_hypertrophy = parse_number(COMP.HYPERTROPHY)) |>
  select(subject = Subject_ID, comp_hypertrophy)

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
    d_myovision_fcsa_I = MyoVision_fCSA_Type_I__Pre_t2 -
      MyoVision_fCSA_Type_I__Pre_t1,
    d_mcsa = mCSA_Pre_t2 - mCSA_Pre_t1,
    d_1rm_legpress = X1RM_Leg_Pre_t2 - X1RM_Leg_Pre_t1,
    d_1rm_ext = `X1RM._Ext_Pre_t2` - `X1RM._Ext_Pre_t1`
  ) |>
  select(subject, all_of(delta_names))

phenotype <- arm |>
  left_join(comp, by = "subject") |>
  left_join(deltas, by = "subject")

cat("Non-missing subjects per trait:\n")
trait_vars <- c("comp_hypertrophy", delta_names)
counts <- vapply(
  trait_vars,
  \(v) sum(!is.na(phenotype[[v]])),
  integer(1L)
)
print(counts)

write_csv(phenotype, here(
  "05_Clustering", "d_integration", "c_data", "phenotype.csv"
))

cat("Done:", nrow(phenotype), "subjects,", ncol(phenotype), "columns\n")
