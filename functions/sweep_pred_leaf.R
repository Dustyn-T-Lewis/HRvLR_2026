# Prediction leaf drivers for the classification and continuous roots. Each leaf
# is one feature level x config x method: a single nested-LOSO cell scored
# against its own permutation null, swept over B. Every leaf stands alone -- its
# workbook holds the observed metric, the per-B permutation p, the null draws,
# the out-of-fold predictions, and the fold selection frequency; its panel shows
# the observed metric against the null band, its selected features beside it.
# No cross-cell comparison lives in any leaf.

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

# Top selected features by fold frequency, or a note when the learner is dense.
selection_panel <- function(selection, fill) {
  top <- selection |> slice_max(.data$freq, n = 8, with_ties = FALSE)
  if (!nrow(top)) {
    return(
      ggplot() +
        annotate("text", 0, 0,
          label = "no recurrent\nfeature selected",
          size = 3, colour = "grey45", fontface = "italic"
        ) +
        theme_void()
    )
  }
  ggplot(top, aes(.data$freq, stats::reorder(.data$feature, .data$freq))) +
    geom_col(fill = fill, width = 0.7) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(x = "fold selection freq.", y = NULL, subtitle = "top features") +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

# Classification leaf: observed AUC against its permutation null with the top
# selected features beside it.
build_class_leaf_panel <- function(res, level, config, method) {
  fill <- LEVEL_FILL[[level]]
  null <- res$null
  obs <- res$summary$estimate[[1]]
  bmax <- max(res$summary$B)
  ptab <- res$summary |>
    filter(.data$B > 0) |>
    mutate(lab = sprintf("B=%d: p=%.3f", .data$B, .data$perm_p))
  lab <- paste(ptab$lab, collapse = "\n")

  p_null <- ggplot(null, aes(.data$auc)) +
    geom_histogram(
      bins = 30, fill = "grey75", colour = "white", linewidth = 0.1
    ) +
    geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey45") +
    geom_vline(xintercept = obs, colour = "#B2182B", linewidth = 0.9) +
    annotate("text", -Inf, Inf,
      label = sprintf("AUC = %.2f\n%s", obs, lab),
      hjust = -0.08, vjust = 1.2, size = 2.5, fontface = "bold",
      colour = "grey20"
    ) +
    labs(
      x = sprintf("permutation-null AUC (B = %d)", bmax), y = "count",
      subtitle = "observed (red) vs null (grey); dotted = chance"
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))

  (p_null | selection_panel(res$selection, fill)) +
    plot_layout(widths = c(1.4, 1)) +
    plot_annotation(
      title = sprintf(
        "%s | %s classify HR/LR -- %s", level, config, MODEL_LABEL[[method]]
      ),
      subtitle = sprintf(
        "%s. Nested LOSO, seeded; permutation null. %s",
        CONFIG_ROLE[[config]], LEAK_NOTE[[level]]
      ),
      theme = leaf_title_theme()
    )
}

# Continuous leaf: observed Q^2 vs null mean per outcome, plus the selected
# features for the outcome with the strongest Q^2.
build_cont_leaf_panel <- function(summ, null, selection, level, config,
                                  method) {
  fill <- LEVEL_FILL[[level]]
  top_b <- summ |> filter(.data$B == max(.data$B))
  disp <- top_b |>
    mutate(
      q2_disp = pmin(pmax(.data$q2, -0.6), 0.9),
      lab = sprintf("%.2f\np=%.3f", .data$q2, .data$perm_p_q2)
    )
  p_q2 <- ggplot(
    disp, aes(stats::reorder(.data$outcome, .data$q2), .data$q2_disp)
  ) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.4) +
    geom_point(aes(y = .data$null_q2_mean),
      shape = 18, size = 2.4, colour = "grey40"
    ) +
    geom_col(fill = fill, width = 0.6, alpha = 0.85) +
    geom_text(aes(y = pmax(.data$q2_disp, 0) + 0.06, label = .data$lab),
      size = 2, lineheight = 0.85, fontface = "bold", colour = "grey20"
    ) +
    scale_y_continuous(limits = c(-0.62, 1)) +
    coord_flip() +
    labs(
      x = NULL, y = expression(LOSO ~ Q^2),
      subtitle = "Q^2 <= 0 = no better than the mean; diamond = null mean"
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))

  best <- summ |>
    filter(.data$B == max(.data$B)) |>
    slice_max(.data$q2, n = 1, with_ties = FALSE) |>
    pull(.data$outcome)
  sel_best <- if ("outcome" %in% names(selection)) {
    filter(selection, .data$outcome == best)
  } else {
    selection[0, , drop = FALSE]
  }

  (p_q2 | selection_panel(sel_best, fill)) +
    plot_layout(widths = c(1.5, 1)) +
    plot_annotation(
      title = sprintf(
        "%s | %s predict adaptation -- %s", level, config, MODEL_LABEL[[method]]
      ),
      subtitle = sprintf(
        "Six deltas, nested LOSO vs permutation null. %s", LEAK_NOTE[[level]]
      ),
      theme = leaf_title_theme()
    )
}

leaf_title_theme <- function() {
  theme(
    plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle = element_text(
      face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
    )
  )
}

LEVEL_FILL <- c(
  proteins = "#2166AC", pathways = "#238B45", modules = "#9E9AC8"
)
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

  dir <- sweep_leaf_dir(root, level, config, method)
  write_sweep_workbook(
    file.path(dir, "c_data", "results.xlsx"),
    list(
      summary = res$summary, null = res$null,
      predictions = res$preds, selection = res$selection
    )
  )
  panel <- build_class_leaf_panel(res, level, config, method)
  save_panel(panel, file.path(dir, "b_reports", "panel"),
    width = 210, height = 95
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

  dir <- sweep_leaf_dir(root, level, config, method)
  write_sweep_workbook(
    file.path(dir, "c_data", "results.xlsx"),
    list(
      summary = summ, null = null, predictions = preds, selection = selection
    )
  )
  panel <- build_cont_leaf_panel(summ, null, selection, level, config, method)
  save_panel(panel, file.path(dir, "b_reports", "panel"),
    width = 220, height = 110
  )
  summ
}
