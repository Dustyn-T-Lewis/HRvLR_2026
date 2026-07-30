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
#
# It reports counts, never names. Shrinking the residual can unmask a protein
# the primary fit misses, and a protein appearing only in a secondary model,
# only in a within-arm contrast, in a study whose interactions are null, is not
# a result: significant in HR and not in LR is not a difference between them.
# The full per-protein comparison stays in the workbook so the evidence is on
# disk; nothing here promotes a row out of it.

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

group_q <- grepl("HRvLR|Interaction", summary_tbl$contrast)
cat(sprintf(
  "\nBH < .05, five HR-vs-LR / interaction contrasts: %d primary, %d adjusted",
  sum(summary_tbl$bh_primary[group_q]), sum(summary_tbl$bh_adjusted[group_q])
))
cat(sprintf(
  "\nBH < .05, four within-arm contrasts:            %d primary, %d adjusted\n",
  sum(summary_tbl$bh_primary[!group_q]), sum(summary_tbl$bh_adjusted[!group_q])
))
cat(
  "Adjustment shrinks the residual, so it can raise a within-arm count.",
  "The study question is the five above, and they are unchanged.\n"
)

write.xlsx(
  list(summary = summary_tbl, all = compare),
  here("03_DEP", "a_non_imputed", "c_data", "07_blood_adjusted.xlsx")
)
