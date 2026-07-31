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

LEVEL_BOOKS <- c(
  proteins = "protein_contrasts.xlsx", pathways = "pathway_contrasts.xlsx",
  modules = "module_contrasts.xlsx"
)

CONTRAST_FAMILY <- c(
  Baseline_HRvLR = "HR - LR", Trained_HRvLR = "HR - LR",
  Acute_HRvLR = "HR - LR",
  Training_Interaction = "interaction", Acute_Interaction = "interaction",
  Training_HR = "within-arm", Training_LR = "within-arm",
  Acute_HR = "within-arm", Acute_LR = "within-arm"
)

# The association half of the story. F05 and F06 answer "can it predict"; F04
# answers "is anything different at all", and the honest summary of that is the
# smallest q each contrast reached, not a count of zeroes.
association_rows <- function() {
  purrr::imap_dfr(LEVEL_BOOKS, function(book, lv) {
    read.xlsx(here("04_Features", lv, "c_data", book), "contrasts") |>
      summarise(
        min_q = min(.data$bh, na.rm = TRUE),
        n_bh = sum(.data$bh < 0.05, na.rm = TRUE),
        .by = "contrast"
      ) |>
      mutate(level = lv)
  }) |>
    mutate(
      level = factor(.data$level, levels = SWEEP_LEVELS),
      family = factor(
        CONTRAST_FAMILY[.data$contrast],
        levels = c("HR - LR", "interaction", "within-arm")
      ),
      contrast = factor(.data$contrast, levels = rev(names(CONTRAST_FAMILY)))
    )
}

association_panel <- function() {
  rows <- association_rows()
  ggplot(rows, aes(.data$min_q, .data$contrast)) +
    geom_vline(
      xintercept = 0.05, linetype = "dashed", colour = "#B2182B",
      linewidth = 0.4
    ) +
    geom_point(aes(shape = .data$family), size = 2.2, fill = "grey92") +
    facet_grid(cols = vars(.data$level)) +
    scale_shape_manual(values = c(21, 22, 23), name = NULL) +
    scale_x_log10(
      limits = c(0.001, 1), breaks = c(0.001, 0.01, 0.05, 0.1, 0.5, 1),
      labels = c(".001", ".01", ".05", ".1", ".5", "1")
    ) +
    labs(
      title = "Association: the nine contrasts at three feature levels",
      subtitle = str_wrap(
        paste(
          "Smallest BH q each contrast reached. Protein and module levels",
          "clear nothing; the pathway points left of .05 are within-arm",
          "contrasts or sets scored on under 15 detected members."
        ), 108
      ),
      x = "smallest BH q in the contrast", y = NULL
    ) +
    FIG_THEME +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      plot.subtitle = element_text(
        face = "italic", size = FIG_SUBTITLE_SIZE - 1, colour = "grey30"
      )
    )
}

KEEPER_ROOTS <- list(
  list(
    root = "F05_classification", kind = "class", metric = "estimate",
    p = "perm_p", null_mean = "null_mean", null_sd = "null_sd",
    label = "Classification: HR vs LR", axis = "AUC",
    note = paste(
      "AUC 0 means no signal, not inverted signal: a model that shrinks to the",
      "intercept predicts the training class proportion, and leave-one-out",
      "shifts that proportion against the held-out subject."
    )
  ),
  list(
    root = "F06_prediction", kind = "cont", metric = "q2",
    p = "perm_p_q2", null_mean = "null_q2_mean", null_sd = "null_q2_sd",
    label = "Prediction: continuous adaptation", axis = "Q2",
    note = "Q2 below 0 is worse than predicting the outcome's own mean."
  )
)

# The null band is an empirical quantile interval, not mean +/- 2 SD.
#
# Q2 is unbounded below, so a permuted fit can reach -300 and one such draw
# inflates the SD past any use: for pathways|T1|svm the null runs -5.4 to 0.54,
# mean + 2 SD puts the upper edge at 1.18, and the empirical 95th percentile is
# 0.258. The observed 0.621 beats all 200 draws. Drawn as an SD band that cell
# sits inside its own null while being coloured as clearing it, which is the
# figure contradicting its own statistic. Quantiles match what perm_p measures.
keeper_null_band <- function(root, row) {
  path <- file.path(
    sweep_leaf_dir(root, row$level, row$config, row$model),
    "c_data", "results.xlsx"
  )
  if (!file.exists(path) || !"null" %in% getSheetNames(path)) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  nulls <- read.xlsx(path, "null")
  if (!is.null(row$outcome) && "outcome" %in% names(nulls)) {
    nulls <- nulls[nulls$outcome == row$outcome, ]
  }
  draws <- nulls[[setdiff(names(nulls), c("outcome", "model"))[1]]]
  if (!length(draws)) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  stats::setNames(
    stats::quantile(draws, c(0.05, 0.95), na.rm = TRUE, names = FALSE),
    c("lo", "hi")
  )
}

# Best is the highest metric, not the lowest p. A reader cherry-picking this
# screen takes the biggest number, so that is the number the null must answer.
keeper_rows <- function(spec) {
  rows <- read.xlsx(
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
      null_lo = NA_real_, null_hi = NA_real_,
      baseline = LEAD_BASELINE[[spec$kind]],
      lead = is_lead(.data$metric, .data$perm_p, spec$kind),
      # The same model and config wins at more than one feature level, so the
      # key has to carry the level or the factor levels collide. The label the
      # reader sees stays clean; the facet already names the level.
      row_key = sprintf("%s|%s|%s", .data$level, .data$model, .data$config),
      row_label = sprintf("%s  (%s)", .data$model, .data$config)
    ) |>
    arrange(.data$level, dplyr::desc(.data$metric))

  band <- vapply(
    seq_len(nrow(rows)), function(i) keeper_null_band(spec$root, rows[i, ]),
    numeric(2)
  )
  rows$null_lo <- band["lo", ]
  rows$null_hi <- band["hi", ]
  rows
}

keeper_panel <- function(spec) {
  rows <- keeper_rows(spec)
  rows$row_key <- factor(rows$row_key, levels = rev(rows$row_key))
  labels <- stats::setNames(rows$row_label, as.character(rows$row_key))

  # Without nulls there is nothing for a cell to clear, and printing "0 of 22
  # clear" would read as a result rather than as an absent test. The panel says
  # which it is.
  has_null <- any(!is.na(rows$null_mean))

  # Hold the view to the points and their nearer band edge.
  span <- range(c(rows$metric, rows$baseline[1]), na.rm = TRUE)
  pad <- max(diff(span), 0.2) * 0.35
  x_limits <- c(span[1] - pad, span[2] + pad)
  n_clipped <- sum(
    rows$null_lo < x_limits[1] | rows$null_hi > x_limits[2],
    na.rm = TRUE
  )

  subtitle <- if (has_null) {
    paste(
      "Best cell per model and feature level against its own permutation null",
      sprintf(
        "(grey, 5th-95th percentile) and the cells it came from. %d of %d",
        sum(rows$lead, na.rm = TRUE), nrow(rows)
      ),
      "clear both the null and the baseline.",
      if (n_clipped) {
        sprintf(
          "%d null band(s) extend past the axis and are clipped, not dropped.",
          n_clipped
        )
      } else {
        ""
      }
    )
  } else {
    paste(
      "Best cell per model and feature level, with the cells it was picked",
      "from. PERMUTATION NULLS NOT YET RUN: no cell can be called a lead, and",
      "the dashed line is the trivial baseline, not a significance threshold."
    )
  }

  subtitle <- paste(subtitle, spec$note)

  ggplot(rows, aes(.data$metric, .data$row_key)) +
    geom_vline(
      xintercept = rows$baseline[1], linetype = "dashed",
      colour = "grey45", linewidth = 0.4
    ) +
    geom_segment(
      aes(x = .data$null_lo, xend = .data$null_hi, yend = .data$row_key),
      colour = "grey72", linewidth = 2.6, lineend = "round"
    ) +
    geom_point(aes(fill = .data$lead), shape = 21, size = 2.6, stroke = 0.4) +
    geom_text(
      aes(x = Inf, label = sprintf("of %d", .data$n_screened)),
      hjust = 1.1, size = 2.1, colour = "grey40"
    ) +
    facet_grid(rows = vars(.data$level), scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = labels) +
    scale_fill_manual(
      values = c(`TRUE` = "#B2182B", `FALSE` = "grey92"),
      labels = c(`TRUE` = "clears null and baseline", `FALSE` = "does not"),
      name = NULL, na.value = "grey92",
      guide = if (has_null) "legend" else "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.16))) +
    # An unpenalised model on a permuted outcome can reach Q2 = -300, and one
    # such band flattens every real point onto a line at zero. Clip rather than
    # drop: coord_cartesian keeps the data and truncates the view, and the
    # subtitle reports how many bands run past the edge.
    coord_cartesian(xlim = x_limits) +
    labs(
      title = spec$label, subtitle = str_wrap(subtitle, 108),
      x = spec$axis, y = NULL
    ) +
    FIG_THEME +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      plot.subtitle = element_text(
        face = "italic", size = FIG_SUBTITLE_SIZE - 1, colour = "grey30"
      )
    )
}
