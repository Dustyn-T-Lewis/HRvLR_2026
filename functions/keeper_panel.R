# The best cell each model family reached, next to the null it was drawn from.
#
# F05 and F06 are screens: 153 classification cells and 792 prediction cells.
# The maximum of a screen estimates nothing, and reporting it alone is how a
# 945-cell sweep becomes a finding. Every row carries the three things a bare
# maximum hides - the permutation null band for that cell, the trivial baseline
# the metric has to beat, and how many cells the winner was picked from.
#
# A point sitting inside its own grey band is noise that happened to come first,
# and the figure should say so at a glance rather than in a footnote.

pacman::p_load(here, dplyr, tidyr, ggplot2, openxlsx, stringr)
source(here("functions", "shared_style.R"))
source(here("functions", "sweep_grid.R"))

KEEPER_ROOTS <- list(
  list(
    root = "F05_classification", kind = "class", metric = "estimate",
    p = "perm_p", null_mean = "null_mean", null_sd = "null_sd",
    label = "Classification: HR vs LR", axis = "AUC"
  ),
  list(
    root = "F06_prediction", kind = "cont", metric = "q2",
    p = "perm_p_q2", null_mean = "null_q2_mean", null_sd = "null_q2_sd",
    label = "Prediction: continuous adaptation", axis = "Q2"
  )
)

# Best is the highest metric, not the lowest p. A reader cherry-picking this
# screen takes the biggest number, so that is the number the null must answer.
keeper_rows <- function(spec) {
  read.xlsx(
    file.path(sweep_root_dir(spec$root), "c_data", "results.xlsx"), "all_cells"
  ) |>
    best_b_per_cell() |>
    mutate(
      metric = .data[[spec$metric]],
      perm_p = .data[[spec$p]],
      null_mean = .data[[spec$null_mean]],
      null_sd = .data[[spec$null_sd]]
    ) |>
    filter(!is.na(.data$metric)) |>
    mutate(n_screened = dplyr::n(), .by = c("level", "model")) |>
    slice_max(.data$metric,
      n = 1, by = c("level", "model"), with_ties = FALSE
    ) |>
    mutate(
      level = factor(.data$level, levels = SWEEP_LEVELS),
      null_lo = .data$null_mean - 2 * .data$null_sd,
      null_hi = .data$null_mean + 2 * .data$null_sd,
      baseline = LEAD_BASELINE[[spec$kind]],
      lead = is_lead(.data$metric, .data$perm_p, spec$kind),
      row_label = sprintf("%s  (%s)", .data$model, .data$config)
    ) |>
    arrange(.data$level, dplyr::desc(.data$metric))
}

keeper_panel <- function(spec) {
  rows <- keeper_rows(spec)
  rows$row_label <- factor(rows$row_label, levels = rev(rows$row_label))
  n_lead <- sum(rows$lead, na.rm = TRUE)
  subtitle <- paste(
    "Best cell per model and feature level against its own permutation null",
    sprintf(
      "(grey, mean +/- 2 SD) and the cells it was picked from. %d of %d clear",
      n_lead, nrow(rows)
    ),
    "both the null and the baseline."
  )

  ggplot(rows, aes(.data$metric, .data$row_label)) +
    geom_vline(
      xintercept = rows$baseline[1], linetype = "dashed",
      colour = "grey45", linewidth = 0.4
    ) +
    geom_segment(
      aes(x = .data$null_lo, xend = .data$null_hi, yend = .data$row_label),
      colour = "grey72", linewidth = 2.6, lineend = "round"
    ) +
    geom_point(aes(fill = .data$lead), shape = 21, size = 2.6, stroke = 0.4) +
    geom_text(
      aes(x = .data$null_hi, label = sprintf("of %d", .data$n_screened)),
      hjust = -0.25, size = 2.1, colour = "grey40"
    ) +
    facet_grid(rows = vars(.data$level), scales = "free_y", space = "free_y") +
    scale_fill_manual(
      values = c(`TRUE` = "#B2182B", `FALSE` = "grey92"),
      labels = c(`TRUE` = "clears null and baseline", `FALSE` = "does not"),
      name = NULL, na.value = "grey92"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    labs(title = spec$label, subtitle = subtitle, x = spec$axis, y = NULL) +
    FIG_THEME +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      plot.subtitle = element_text(
        face = "italic", size = FIG_SUBTITLE_SIZE - 1, colour = "grey30"
      )
    )
}
