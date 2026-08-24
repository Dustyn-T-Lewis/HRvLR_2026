# Panel C: what the sweep produced against what noise produces. Every count
# falls below its expectation, which is the result: the cells that cleared BH
# are fewer than a 180-cell sweep yields with nothing in it.
if (!exists("sweep_calibration")) {
  source(here::here("04_Figures", "F03_modules", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, tidyr, dplyr, forcats)

build_calibration <- function(sweep_calibration, tag = "C") {
  d <- sweep_calibration |>
    filter(.data$level == "all") |>
    dplyr::select(
      `cells with a hit` = observed_hit_cells,
      `total hits` = observed_hits,
      `cells at p < 0.05` = cells_p_below_05
    ) |>
    pivot_longer(everything(), names_to = "stat", values_to = "observed") |>
    left_join(
      sweep_calibration |>
        filter(.data$level == "all") |>
        dplyr::select(
          `cells with a hit` = expected_hit_cells,
          `total hits` = expected_hits,
          `cells at p < 0.05` = n_cells
        ) |>
        mutate(`cells at p < 0.05` = 0.05 * .data$`cells at p < 0.05`) |>
        pivot_longer(everything(), names_to = "stat", values_to = "expected"),
      by = "stat"
    ) |>
    pivot_longer(c("observed", "expected"), names_to = "source") |>
    mutate(stat = fct_inorder(.data$stat))

  p <- ggplot(d, aes(.data$value, .data$stat, fill = .data$source)) +
    geom_col(position = position_dodge(0.7), width = 0.6) +
    geom_text(aes(label = round(.data$value, 1)),
      position = position_dodge(0.7), hjust = -0.25, size = 2.2
    ) +
    scale_fill_manual(
      values = c(observed = "#2166AC", expected = "grey70"),
      breaks = c("observed", "expected"), name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      title = "The sweep is quieter than chance", tag = tag,
      subtitle = paste(
        "180 cells. Expectation from shuffling the phenotype,",
        "999 times per cell"
      ),
      x = NULL, y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      legend.position = "top",
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = sweep_calibration)
}

c_calibration <- build_calibration(sweep_calibration)
save_panel(
  c_calibration$plot, file.path(F03_RPT, "panels", "panel_c_calibration"),
  150, 80
)
F03_PANELS[["calibration"]] <- c_calibration$plot
F03_AUDIT[["sweep_calibration"]] <- c_calibration$audit
