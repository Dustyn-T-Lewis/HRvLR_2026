# Do the prediction leads survive with no imputation at all?
#
# Every F06 cell is fitted on the missForest matrix, imputed once across all 45
# samples, so a held-out subject's filled values were informed by the subjects
# training the model that predicts them. loso_refit already showed that the
# module leads collapse when their other leak, cohort-defined eigengenes, is
# closed. The pathway and protein leads have never been tested against this one.
#
# The test refits each lead on proteins observed in every sample, with nothing
# imputed. A lead that holds was not built on filled values. Modules cannot be
# tested this way and are not attempted: an eigengene needs the WGCNA network,
# and that network is the separate leak loso_refit measures.

pacman::p_load(here, dplyr, tibble, openxlsx)
source(here("functions", "pred_complete_case.R"))
source(here("functions", "shared_prediction.R"))
source(here("functions", "sweep_grid.R"))

B <- 200L

cells <- read.xlsx(
  here("05_Figures", "F06_prediction", "c_data", "results.xlsx"), "all_cells"
) |>
  best_b_per_cell() |>
  filter(
    is_lead(.data$q2, .data$perm_p_q2, "cont"),
    .data$level %in% c("pathways", "proteins")
  )

message(sprintf("rescoring %d non-module leads without imputation", nrow(cells)))

bundle <- pred_load_complete_case()
key <- c(pathways = "singscore", proteins = "proteins")

rescored <- bind_rows(lapply(seq_len(nrow(cells)), function(i) {
  cell <- cells[i, ]
  x <- pred_contrast_matrix(
    bundle$feature_sets[[key[[cell$level]]]], bundle$meta, cell$config
  )
  al <- align_xy(x, pred_outcome(bundle, cell$outcome))
  r <- sweep_cont_cell(al$x, al$y, cell$model, cell$outcome, b_grid = c(B))
  s <- r$summary[r$summary$B == B, ]
  tibble(
    level = cell$level, config = cell$config, outcome = cell$outcome,
    model = cell$model, n_features = ncol(al$x),
    q2_imputed = cell$q2, p_imputed = cell$perm_p_q2,
    q2_complete = s$q2, p_complete = s$perm_p_q2
  )
})) |>
  mutate(
    drop = .data$q2_imputed - .data$q2_complete,
    survives = .data$q2_complete > 0 & .data$p_complete < 0.05
  ) |>
  arrange(dplyr::desc(.data$q2_imputed))

print(as.data.frame(rescored), digits = 3)

cat(sprintf(
  "\n%d of %d leads survive with no imputation. Median Q2 drop %.3f.\n",
  sum(rescored$survives), nrow(rescored), median(rescored$drop)
))

write.xlsx(
  list(rescored = rescored),
  here("05_Figures", "F06_prediction", "c_data", "imputation_sensitivity.xlsx")
)
