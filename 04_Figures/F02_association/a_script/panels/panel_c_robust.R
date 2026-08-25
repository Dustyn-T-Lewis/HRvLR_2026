# Panel C: refit the survivor fifteen times, dropping one subject each time.
# A result that holds in six of fifteen folds is a result about whichever
# subjects were kept. The rank correlation in the subtitle is the second
# check: a linear fit that a Spearman test cannot see is being carried by the
# extremes of the scale rather than by the ordering.
if (!exists("F02_RPT")) {
  source(here::here("04_Figures", "F02_association", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats, dplyr, openxlsx)

build_robust <- function(tag = "C") {
  book <- here::here("03_Features", "c_data", "04_hit_robustness.xlsx")
  rb <- as_tibble(read.xlsx(book, "robustness")) |> slice_min(bh_full, n = 1)
  lo <- as_tibble(read.xlsx(book, "loso_detail")) |>
    filter(.data$feature == rb$feature[1], .data$window == rb$window[1]) |>
    mutate(
      dropped = fct_reorder(.data$dropped, .data$bh),
      holds = .data$bh < 0.05
    )

  p <- ggplot(lo, aes(.data$bh, .data$dropped, colour = .data$holds)) +
    geom_vline(
      xintercept = 0.05, linetype = "dashed", colour = "grey50",
      linewidth = 0.4
    ) +
    geom_segment(aes(x = 0, xend = .data$bh, yend = .data$dropped),
      linewidth = 0.4
    ) +
    geom_point(size = 2.2) +
    geom_point(
      data = rb, aes(x = .data$bh_full, y = Inf), inherit.aes = FALSE,
      colour = "grey30", shape = 18, size = 3
    ) +
    scale_colour_manual(
      values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"),
      labels = c(`TRUE` = "still < 0.05", `FALSE` = "lost"), name = NULL
    ) +
    labs(
      title = "The strongest hit, refit without each subject", tag = tag,
      subtitle = sprintf(
        "Holds in %d of %d folds. Spearman rho = %.2f, p = %.2f",
        rb$loso_folds_below_alpha[1], rb$loso_folds[1],
        rb$spearman_rho[1], rb$spearman_p[1]
      ),
      x = "BH-adjusted p with that subject dropped", y = "Subject dropped"
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      axis.text.y = element_text(size = 6),
      legend.position = "top",
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = lo |> dplyr::select(-holds), summary = rb)
}

c_robust <- build_robust()
save_panel(
  c_robust$plot, file.path(F02_RPT, "panels", "panel_c_robust"), 110, 105
)
F02_PANELS[["robust"]] <- c_robust$plot
F02_AUDIT[["hit_loso"]] <- c_robust$audit
F02_AUDIT[["hit_robustness"]] <- c_robust$summary
