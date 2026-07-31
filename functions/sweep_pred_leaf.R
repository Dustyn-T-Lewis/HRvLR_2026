# Prediction leaf drivers for the classification and continuous roots. Each leaf
# is one feature level x config x method: a single nested-LOSO cell scored
# against its own permutation null, swept over B. Every leaf stands alone -- the
# store holds its observed metric, the per-B permutation p, the null draws, the
# out-of-fold predictions, and the fold selection frequency, one row set per
# phenotype. No cross-cell comparison lives in any leaf.

pacman::p_load(here, dplyr, tidyr, ggplot2, patchwork, scales)
source(here("functions", "shared_style.R"))
source(here("functions", "sweep_grid.R"))
source(here("functions", "pred_features.R"))
source(here("functions", "shared_prediction.R"))

MODEL_LABEL <- c(
  enet = "elastic net", lasso = "lasso", ridge = "ridge", spls = "sPLS",
  pam = "PAM", rf = "random forest", svm = "linear SVM", plain = "plain"
)

leaf_x <- function(bundle, level, config) {
  key <- SWEEP_LEVEL_KEY[[level]]
  pred_contrast_matrix(bundle$feature_sets[[key]], bundle$meta, config)
}

leaf_title_theme <- function() {
  theme(
    plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle = element_text(
      face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
    )
  )
}

LEAK_NOTE <- c(
  proteins = "proteins cohort-imputed (optimistic).",
  pathways = "pathways leakage-free.",
  modules = "module eigengenes cohort-relative (optimistic)."
)

run_class_sweep_leaf <- function(bundle, level, config, method, b_grid = B_GRID,
                                 root = "F05_classification") {
  x0 <- leaf_x(bundle, level, config)
  al <- align_xy(x0, pred_outcome(bundle, "group"))
  res <- sweep_class_cell(al$x, al$y, method, b_grid = b_grid)
  res$summary <- cbind(level = level, config = config, res$summary)

  write_sweep_cell(
    root, level, config, "HR_LR", method,
    list(
      summary = res$summary, null = res$null,
      predictions = res$preds, selection = res$selection
    ),
    fingerprint = sweep_fingerprint(bundle, b_grid)
  )
  res$summary
}

run_cont_sweep_leaf <- function(bundle, level, config, method, b_grid = B_GRID,
                                outcomes = ADAPT_OUTCOMES,
                                root = "F06_prediction") {
  x0 <- leaf_x(bundle, level, config)
  cells <- lapply(outcomes, function(oc) {
    al <- align_xy(x0, pred_outcome(bundle, oc))
    r <- sweep_cont_cell(al$x, al$y, method, oc, b_grid = b_grid)
    if (nrow(r$selection)) {
      r$selection <- cbind(outcome = oc, r$selection)
    }
    r
  })
  summ <- bind_rows(lapply(cells, `[[`, "summary")) |>
    mutate(level = level, config = config, .before = 1)
  null <- bind_rows(lapply(cells, `[[`, "null"))
  preds <- bind_rows(lapply(cells, `[[`, "preds"))
  selection <- bind_rows(lapply(cells, `[[`, "selection"))

  # One cell per outcome, sliced here rather than by a later pass over the
  # written file.
  fingerprint <- sweep_fingerprint(bundle, b_grid)
  for (oc in outcomes) {
    write_sweep_cell(
      root, level, config, oc, method,
      list(
        summary = filter_to_outcome(summ, oc),
        null = filter_to_outcome(null, oc),
        predictions = filter_to_outcome(preds, oc),
        selection = filter_to_outcome(selection, oc)
      ),
      fingerprint = fingerprint
    )
  }
  summ
}
