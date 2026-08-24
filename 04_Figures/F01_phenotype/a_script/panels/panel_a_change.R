# Panel A: did each outcome change over training at all? On a common effect-size
# scale because the raw units span two orders of magnitude and the comparison
# across outcomes is the whole point of the panel.
if (!exists("change_summary")) {
  source(here::here("04_Figures", "F01_phenotype", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats)

build_change <- function(change_summary, tag = "A") {
  d <- change_summary |>
    filter(.data$trait != "comp_hypertrophy") |>
    mutate(
      label = TRAIT_LABELS[.data$trait],
      moved = .data$ci_lo_d > 0 | .data$ci_hi_d < 0,
      label = factor(.data$label, levels = TRAIT_ORDER)
    )

  p <- ggplot(d, aes(.data$mean_d, .data$label, colour = .data$moved)) +
    geom_vline(
      xintercept = 0, linetype = "dashed", colour = "grey50",
      linewidth = 0.4
    ) +
    geom_errorbar(aes(xmin = .data$ci_lo_d, xmax = .data$ci_hi_d),
      orientation = "y", width = 0.18, linewidth = 0.6
    ) +
    geom_point(size = 2.6) +
    scale_colour_manual(
      values = c(`TRUE` = "#2166AC", `FALSE` = "grey60"), guide = "none"
    ) +
    labs(
      title = "Which outcomes actually changed", tag = tag,
      subtitle = "Mean change over training, in SD units, with 95% CI",
      x = "Change (SD units)", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = d |> select(
    trait, n, mean, ci_lo, ci_hi,
    mean_d, ci_lo_d, ci_hi_d, p, moved
  ))
}

a_change <- build_change(change_summary)
save_panel(
  a_change$plot, file.path(F01_RPT, "panels", "panel_a_change"),
  135, 85
)
F01_PANELS[["change"]] <- a_change$plot
F01_AUDIT[["change_over_training"]] <- a_change$audit
