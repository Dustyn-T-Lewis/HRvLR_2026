# Can the proteome classify HR vs LR out of sample? Three feature levels
# (protein abundance, singscore, module eigengene) x three timepoint
# windows (T1 baseline, training = T2-T1, acute = T3-T2) x one outcome
# (HR/LR), elastic net via glmnet::cv.glmnet, leave-one-subject-out,
# subject-label permutation null. Plus one plain-glm sensitivity row for
# modules only, where p < n holds (12 modules < 16 subjects) and an
# unpenalized fit is defined.
#
# 12 cells total, against the deleted F05_classification screen's 153 --
# scoped to what this repo's own prior build notes already concluded:
# complex learners don't beat regularized linear models at n=16, so one
# learner, not a menu; T2/T3 raw snapshots and the weakest configs (total,
# trajectory) are dropped, keeping the three windows that carry a real,
# distinct question.
#
# No BH across the 12 cells: this repo's own README already makes the
# argument for its now-deleted screens, and it applies here too -- a
# screen this size and this correlated (shared subjects, overlapping
# feature spaces) has no defensible multiple-comparison correction. Raw
# permutation p plus the exact cell count is the honest alternative. A
# lead needs perm_p < .05 AND to beat the trivial baseline (AUC > 0.5) --
# a permutation null can reach its p floor while still predicting worse
# than chance.

pacman::p_load(here, dplyr, purrr, tibble, openxlsx)
source(here("functions", "pred_features.R"))
source(here("functions", "shared_prediction.R"))

OUT_DIR <- here("03_Analysis", "categorical", "F04_classification", "c_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LEVELS <- c("proteins", "singscore", "eigengenes")
CONFIGS <- c("T1", "training", "acute")

bundle <- pred_load("categorical")
y_all <- pred_outcome(bundle, "group")

run_cell <- function(level, config, model) {
  feat <- pred_contrast_matrix(
    bundle$feature_sets[[level]], bundle$meta, config
  )
  al <- align_xy(feat, y_all)
  res <- suppressWarnings(run_class_cell(al$x, al$y, model))
  res$summary |>
    mutate(
      level = level, config = config, .before = 1,
      lead = perm_p < 0.05 & estimate > 0.5
    )
}

grid <- tidyr::expand_grid(level = LEVELS, config = CONFIGS)
cat(sprintf("Running %d classification cells (enet)...\n", nrow(grid)))

enet_results <- purrr::pmap_dfr(grid, function(level, config) {
  run_cell(level, config, "enet")
})

cat("Running 3 plain-glm sensitivity cells (modules only)...\n")
plain_results <- purrr::map_dfr(CONFIGS, function(config) {
  run_cell("eigengenes", config, "plain") |>
    mutate(model = "plain (sensitivity)")
})

results <- bind_rows(enet_results, plain_results) |>
  select(level, config, model, n, estimate, ci_lo, ci_hi, perm_p, lead)

print(as.data.frame(results), digits = 3)
cat(sprintf(
  "\n%d of %d primary (enet) cells are leads (perm_p < .05 and AUC > 0.5)\n",
  sum(enet_results$perm_p < 0.05 & enet_results$estimate > 0.5), nrow(grid)
))

write.xlsx(
  list(cells = results),
  file.path(OUT_DIR, "F04_classification_source_data.xlsx")
)
