# Specification-curve view of a prediction root (Simonsohn 2015). Every cell in
# the screen is one point, sorted by its observed metric, drawn against its own
# permutation-null band; a strip below reads off the feature level, config, and
# method behind each point. It shows the whole distribution and never ranks one
# cell over another as a result -- a review aid for choosing a method, not a
# claim.

pacman::p_load(here, dplyr, tidyr, ggplot2, patchwork, forcats)
source(here("functions", "shared_style.R"))
source(here("functions", "sweep_grid.R"))

SPEC_LEVEL_COLORS <- c(
  proteins = "#2166AC", pathways = "#238B45", modules = "#9E9AC8"
)

# The read-off axes below the curve: one facet per specification dimension,
# each cell's choice marked at its rank. Shared by all three screens so the
# curves stay comparable.
spec_strip <- function(d, spec_cols) {
  strip <- d |>
    dplyr::select(.data$rank, dplyr::all_of(spec_cols)) |>
    pivot_longer(all_of(spec_cols), names_to = "axis", values_to = "spec") |>
    mutate(axis = factor(axis, levels = spec_cols))

  ggplot(strip, aes(.data$rank, fct_rev(.data$spec))) +
    geom_point(aes(colour = .data$axis), size = 0.9) +
    scale_colour_manual(
      values = stats::setNames(
        c("#2166AC", "#E69F00", "#238B45", "#6A51A3")[seq_along(spec_cols)],
        spec_cols
      ),
      guide = "none"
    ) +
    facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "specification (cells sorted by metric)", y = NULL) +
    FIG_THEME +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 7),
      panel.spacing = unit(1, "pt")
    )
}

# cells: roll-up all_cells rows at one B. value/mean/sd/p name the metric and
# its null; chance is the no-skill reference; spec_cols are the read-off axes.
# y_lim clamps the display so a few catastrophic negatives (a plain model
# overfitting the module space to Q^2 of -1e6) cannot squash the informative
# range; clamped points and bands pile at the floor via scales::squish.
spec_curve <- function(cells, value, mean, sd, p, chance, metric_label,
                       spec_cols, title, subtitle, y_lim = NULL) {
  d <- cells |>
    mutate(
      value = .data[[value]], null_mean = .data[[mean]],
      null_sd = .data[[sd]], perm_p = .data[[p]],
      lo = .data[[mean]] - 1.96 * .data[[sd]],
      hi = .data[[mean]] + 1.96 * .data[[sd]],
      hit = is_lead_at(.data[[value]], .data[[p]], chance)
    ) |>
    arrange(.data$value) |>
    mutate(rank = dplyr::row_number())

  y_scale <- if (is.null(y_lim)) {
    scale_y_continuous()
  } else {
    scale_y_continuous(limits = y_lim, oob = scales::squish)
  }

  p_curve <- ggplot(d, aes(.data$rank, .data$value)) +
    geom_linerange(aes(ymin = .data$lo, ymax = .data$hi),
      colour = "grey80", linewidth = 0.5
    ) +
    geom_point(aes(y = .data$null_mean), colour = "grey55", size = 0.35) +
    geom_hline(yintercept = chance, linetype = "dashed", colour = "grey40") +
    geom_point(aes(colour = level, shape = .data$hit), size = 1.4) +
    scale_colour_manual(values = SPEC_LEVEL_COLORS, name = "level") +
    scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 17), guide = "none"
    ) +
    y_scale +
    labs(
      x = NULL, y = metric_label, title = title, subtitle = subtitle
    ) +
    FIG_THEME +
    theme(
      legend.position = c(0.14, 0.8),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE),
      axis.text.x = element_blank(), axis.ticks.x = element_blank()
    )

  p_curve / spec_strip(d, spec_cols) +
    plot_layout(heights = c(1.5, 1.4)) +
    plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
}

spec_curve_class <- function(cells, title, subtitle) {
  spec_curve(
    best_b_per_cell(cells),
    value = "estimate", mean = "null_mean", sd = "null_sd", p = "perm_p",
    chance = 0.5, metric_label = "LOSO AUC",
    spec_cols = c("level", "config", "model"),
    title = title, subtitle = subtitle
  )
}

spec_curve_cont <- function(cells, title, subtitle) {
  spec_curve(
    best_b_per_cell(cells),
    value = "q2", mean = "null_q2_mean", sd = "null_q2_sd", p = "perm_p_q2",
    chance = 0, metric_label = expression(LOSO ~ Q^2),
    spec_cols = c("level", "config", "model", "outcome"),
    title = title, subtitle = subtitle, y_lim = c(-1, 0.7)
  )
}

# Association has no permutation null, so this curve draws no null band and
# marks no leads. It also cannot rank cells by their metric: limma and lm carry
# a moderated t, wilcoxon carries a group effect, and the two are different
# scales. The nominal hit rate is the one quantity every cell reports on the
# same scale, and 5% is where it lands when nothing is there.
render_root_speccurve <- function(root) {
  cells <- openxlsx::read.xlsx(
    file.path(sweep_root_dir(root), "c_data", "results.xlsx"), "all_cells"
  )
  levels_txt <- paste(sort(unique(cells$level)), collapse = ", ")
  is_class <- root == "F05_classification"
  panel <- if (is_class) {
    spec_curve_class(
      cells, "Classification screen -- specification curve",
      sprintf(
        "HR/LR by nested LOSO. Each cell vs its null band. Levels: %s.",
        levels_txt
      )
    )
  } else {
    spec_curve_cont(
      cells, "Continuous-prediction screen -- specification curve",
      sprintf(
        paste(
          "Six deltas by nested LOSO. Each cell vs its null band; Q^2 clamped",
          "at -1. Chance = 0. Levels: %s."
        ),
        levels_txt
      )
    )
  }
  dir.create(file.path(sweep_root_dir(root), "b_reports"),
    recursive = TRUE, showWarnings = FALSE
  )
  save_panel(
    panel, file.path(sweep_root_dir(root), "b_reports", "specification_curve"),
    width = 300, height = if (is_class) 200 else 230
  )
  invisible(panel)
}
