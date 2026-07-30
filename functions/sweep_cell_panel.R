# One panel per split cell. Left carries the distribution the cell's number sits
# in, right carries that phenotype's own named drivers, and the footer states
# the metric, the raw permutation p, how many cells were screened, and whether
# the level is leakage-free. No BH q reaches a panel.

pacman::p_load(here, dplyr, ggplot2, patchwork, ggrepel)
source(here("functions", "shared_style.R"))
source(here("functions", "sweep_drivers.R"))
source(here("functions", "sweep_pred_leaf.R"))

SCREEN_SIZE <- c(F05_classification = 153L, F06_prediction = 792L)

# `clean_pathway_name` arrives via sweep_drivers.R -> shared_pathway_utils.R.
# Proteins and modules already name themselves.
cell_feature_label <- function(feature, level) {
  if (level == "pathways") clean_pathway_name(feature) else feature
}

cell_footer <- function(metric_lab, metric, p, root, level) {
  sprintf(
    "%s = %.2f  ·  p = %s  ·  1 of %d cells screened  ·  %s",
    metric_lab, metric,
    if (is.na(p)) "not estimated" else sprintf("%.3f", p),
    SCREEN_SIZE[[root]], LEAK_NOTE[[level]]
  )
}

# Ridge retains every coefficient, so its fold frequencies describe the whole
# feature space rather than a signature.
drivers_or_empty <- function(selection, level, method, xlab, title, subtitle,
                             signed = FALSE, n = 10) {
  dense <- method %in% DENSE_MODELS || is.null(selection) || !nrow(selection)
  d <- if (dense) {
    selection[0, , drop = FALSE]
  } else {
    selection |>
      slice_max(abs(.data$score), n = n, with_ties = FALSE)
  }
  driver_bars(d, level, xlab, title, subtitle, signed = signed)
}

# Wilcoxon is rank-based and carries no t, so a cell scores on whatever its
# method produced: the moderated t where it exists, the group difference
# otherwise.
null_panel <- function(null_values, observed, chance, xlab, subtitle) {
  ggplot(data.frame(v = null_values), aes(.data$v)) +
    geom_histogram(
      bins = 30, fill = "grey75", colour = "white", linewidth = 0.1
    ) +
    geom_vline(xintercept = chance, linetype = "dotted", colour = "grey45") +
    geom_vline(xintercept = observed, colour = "#B2182B", linewidth = 0.9) +
    labs(x = xlab, y = "count", subtitle = subtitle) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

build_class_cell_panel <- function(sheets, level, config, method) {
  # sheets$summary is pre-filtered to one phenotype by sweep_split.R, so the
  # max-B row is a single row and `[[1]]` is safe here and below.
  s <- sheets$summary |> filter(.data$B == max(.data$B))
  obs <- s$estimate[[1]]

  p_null <- null_panel(
    sheets$null$auc, obs, 0.5,
    sprintf("permutation-null AUC (B = %d)", s$B[[1]]),
    "observed (red) vs null (grey); dotted = chance"
  )
  bars <- drivers_or_empty(
    sheets$selection |> mutate(score = .data$freq), level, method,
    "fold selection freq.", NULL, "drivers of the HR/LR assignment"
  )

  (p_null | bars) +
    plot_layout(widths = c(1.4, 1)) +
    plot_annotation(
      title = sprintf(
        "%s · %s · HR_LR classify -- %s", level, config, MODEL_LABEL[[method]]
      ),
      subtitle = cell_footer(
        "AUC", obs, s$perm_p[[1]], "F05_classification", level
      ),
      theme = leaf_title_theme()
    )
}

build_cont_cell_panel <- function(sheets, level, config, phenotype, method) {
  # sheets$summary is pre-filtered to one phenotype by sweep_split.R, so the
  # max-B row is a single row and `[[1]]` is safe here and below.
  s <- sheets$summary |> filter(.data$B == max(.data$B))
  obs <- s$q2[[1]]

  p_null <- null_panel(
    sheets$null$q2, obs, 0,
    sprintf("permutation-null Q2 (B = %d)", s$B[[1]]),
    "observed (red) vs null (grey); dotted = predicting the mean"
  )
  bars <- drivers_or_empty(
    sheets$selection |> mutate(score = .data$freq), level, method,
    "fold selection freq.", NULL, sprintf("drivers of %s", phenotype)
  )

  (p_null | bars) +
    plot_layout(widths = c(1.4, 1)) +
    plot_annotation(
      title = sprintf(
        "%s · %s · %s predict -- %s", level, config, phenotype,
        MODEL_LABEL[[method]]
      ),
      subtitle = cell_footer(
        "Q2", obs, s$perm_p_q2[[1]], "F06_prediction", level
      ),
      theme = leaf_title_theme()
    )
}
