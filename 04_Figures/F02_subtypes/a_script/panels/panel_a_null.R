# Panel A: the observed fit gain against what the same procedure produces on
# data with no clusters in it. The null is the panel - without it a BIC gain is
# unreadable at this sample size.
if (!exists("null_draws")) {
  source(here::here("04_Figures", "F02_subtypes", "a_script", "setup.R"))
}
pacman::p_load(ggplot2)

build_null <- function(null_draws, cluster_cells, tag = "A") {
  nd <- null_draws |>
    mutate(
      space_lab = SPACE_LABELS[.data$space],
      cell = paste0(.data$n_pc, " PC")
    )
  cc <- cluster_cells |>
    mutate(
      space_lab = SPACE_LABELS[.data$space],
      cell = paste0(.data$n_pc, " PC")
    )

  p <- ggplot(nd, aes(.data$cell, .data$null_gain)) +
    geom_violin(
      fill = "grey85", colour = "grey60", linewidth = 0.3,
      scale = "width"
    ) +
    geom_point(
      data = cc, aes(.data$cell, .data$bic_gain),
      colour = "#B2182B", size = 2.8
    ) +
    geom_text(
      data = cc, aes(.data$cell, .data$bic_gain,
        label = sprintf("p = %.2f", .data$p_empirical)
      ),
      colour = "#B2182B", size = 2.4, vjust = -1.1
    ) +
    facet_wrap(~space_lab) +
    scale_y_sqrt(breaks = c(0, 10, 50, 100, 200, 400)) +
    labs(
      title = "Observed clustering against a no-cluster null", tag = tag,
      subtitle = paste(
        "Grey = 999 draws from one Gaussian with the observed covariance;",
        "red = observed"
      ),
      x = "Principal components retained",
      y = "BIC gain over one component (sqrt scale)"
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40")
    )

  list(plot = p, audit = cc |> dplyr::select(
    space, n_pc, best_g, bic_gain,
    null_gain_median, null_gain_q95, p_empirical
  ))
}

a_null <- build_null(null_draws, cluster_cells)
save_panel(a_null$plot, file.path(F02_RPT, "panels", "panel_a_null"), 160, 95)
F02_PANELS[["null"]] <- a_null$plot
F02_AUDIT[["cluster_vs_null"]] <- a_null$audit
