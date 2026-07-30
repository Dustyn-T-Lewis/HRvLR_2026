# Sensitivity: the nine contrasts refitted with the per-sample blood index as a
# covariate.
#
# The primary fit is unchanged and stays the headline. This one exists because
# the T3 blood rise is not equal in the two arms - arm x T3 gives b = -1.21,
# p = 0.032, and p = 0.017 by subject-label permutation - so a
# difference-in-differences does not cancel it. Whatever the acute contrasts
# report, they report it against a residual blood difference of roughly 1.6 log2
# units, and this shows how much of that report survives adjustment.
#
# Adjusting is not free: the blood index is built from haemoglobin, and any
# acute biology that moves with perfusion is partly absorbed. That is why it is
# a sensitivity analysis and not the primary.

pacman::p_load(here, dplyr, tibble, openxlsx)
source(here("functions", "feature_contrasts.R"))
source(here("functions", "feature_levels.R"))
source(here("functions", "blood_index_model.R"))

mat <- protein_matrix()
blood <- blood_index_data()
covar <- setNames(blood$blood_index, blood$Col_ID)

primary <- fit_feature_contrasts(mat)
adjusted <- fit_feature_contrasts(mat, adjust = covar)

compare <- primary |>
  select("contrast", "feature", lfc_primary = "logFC", bh_primary = "bh") |>
  inner_join(
    select(adjusted, "contrast", "feature",
      lfc_adjusted = "logFC", bh_adjusted = "bh"
    ),
    by = c("contrast", "feature")
  )

summary_tbl <- compare |>
  group_by(.data$contrast) |>
  summarise(
    n = dplyr::n(),
    bh_primary = sum(.data$bh_primary < 0.05, na.rm = TRUE),
    bh_adjusted = sum(.data$bh_adjusted < 0.05, na.rm = TRUE),
    min_bh_adjusted = min(.data$bh_adjusted, na.rm = TRUE),
    median_lfc_shrink = median(
      abs(.data$lfc_adjusted) - abs(.data$lfc_primary),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(as.data.frame(summary_tbl), digits = 3)

survivors <- compare |>
  filter(.data$bh_primary < 0.05 | .data$bh_adjusted < 0.05) |>
  arrange(.data$bh_adjusted)

if (nrow(survivors)) {
  cat("\nProteins clearing BH in either fit:\n")
  print(as.data.frame(survivors), digits = 3)
} else {
  cat("\nNo protein clears BH in either fit.\n")
}

write.xlsx(
  list(summary = summary_tbl, survivors = survivors, all = compare),
  here("03_DEP", "a_non_imputed", "c_data", "07_blood_adjusted.xlsx")
)
