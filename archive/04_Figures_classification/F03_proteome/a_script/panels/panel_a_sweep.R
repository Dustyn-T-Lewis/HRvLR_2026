# Panel A: how close each split-contrast cell came. Plotting the smallest
# adjusted p rather than a survivor count, because 22 of the 24 cells have a
# count of zero and a bar chart of zeros says nothing about how far from the
# line they fell. One row per partition, not per label: the three fibre-area
# measures cut the cohort identically and are one test, not three.
if (!exists("sweep_summary")) {
  source(here::here("04_Figures", "F03_proteome", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats, dplyr)

build_sweep <- function(sweep_summary, tag = "A") {
  d <- sweep_summary |>
    distinct(.data$partition, .data$contrast, .keep_all = TRUE) |>
    left_join(partition_names(sweep_summary), by = "partition") |>
    mutate(
      label_lab = fct_rev(fct_inorder(.data$split)),
      origin = ifelse(.data$internal, "in composite", "independent"),
      neglog = -log10(.data$min_bh)
    )

  p <- ggplot(d, aes(.data$neglog, .data$label_lab,
    colour = .data$contrast,
    shape = .data$origin
  )) +
    geom_vline(
      xintercept = -log10(0.05), linetype = "dashed",
      colour = "grey50", linewidth = 0.4
    ) +
    geom_point(size = 2.8, position = position_dodge(0.5)) +
    scale_colour_manual(values = unname(CONTRAST_COLORS)[1:3], name = NULL) +
    scale_shape_manual(values = c(
      `in composite` = 17,
      `independent` = 16
    ), name = NULL) +
    labs(
      title = "How close each split came", tag = tag,
      subtitle = sprintf(
        paste(
          "Smallest BH-adjusted p per split and contrast; dashed = 0.05.",
          "%d of %d cells cleared it, against %.1f expected by chance"
        ),
        sweep_calibration$observed_hit_cells[1],
        sweep_calibration$n_cells[1],
        sweep_calibration$expected_hit_cells[1]
      ),
      x = expression(-log[10] * "(smallest adjusted p)"), y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      legend.position = "top", legend.box = "vertical",
      legend.spacing.y = unit(1, "pt"),
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = d |> dplyr::select(
    partition, split, contrast, internal, n_tested,
    n_nominal, n_bh, min_bh
  ))
}

a_sweep <- build_sweep(sweep_summary)
save_panel(
  a_sweep$plot, file.path(F03_RPT, "panels", "panel_a_sweep"),
  150, 110
)
F03_PANELS[["sweep"]] <- a_sweep$plot
F03_AUDIT[["sweep_summary"]] <- a_sweep$audit
