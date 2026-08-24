# Panel A: which adaptations actually happened. On a common effect-size scale
# because the raw units span orders of magnitude, and because an outcome whose
# change interval covers zero cannot anchor an association no matter what the
# proteome does.
if (!exists("pheno")) {
  source(here::here("04_Figures", "F01_phenotype", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats, broom, purrr, dplyr)

build_change <- function(pheno, tag = "A") {
  d <- map_dfr(CHANGE_TRAITS, function(v) {
    broom::tidy(stats::t.test(pheno[[v]])) |>
      transmute(
        trait = v, n = sum(!is.na(pheno[[v]])),
        sd = stats::sd(pheno[[v]], na.rm = TRUE),
        mean_d = estimate / sd, ci_lo_d = conf.low / sd,
        ci_hi_d = conf.high / sd, p = p.value
      )
  }) |>
    mutate(
      label = fct_reorder(TRAIT_LABELS[.data$trait], .data$mean_d),
      moved = .data$ci_lo_d > 0 | .data$ci_hi_d < 0
    )

  p <- ggplot(d, aes(.data$mean_d, .data$label, colour = .data$moved)) +
    geom_vline(
      xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4
    ) +
    geom_errorbar(aes(xmin = .data$ci_lo_d, xmax = .data$ci_hi_d),
      orientation = "y", width = 0.18, linewidth = 0.6
    ) +
    geom_point(size = 2.6) +
    scale_colour_manual(
      values = c(`TRUE` = "#2166AC", `FALSE` = "grey60"), guide = "none"
    ) +
    labs(
      title = "Which adaptations happened", tag = tag,
      subtitle = "Mean change over training, in SD units, with 95% CI",
      x = "Change (SD units)", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = d |> dplyr::select(-label))
}

a_change <- build_change(pheno)
save_panel(
  a_change$plot, file.path(F01_RPT, "panels", "panel_a_change"),
  135, 90
)
F01_PANELS[["change"]] <- a_change$plot
F01_AUDIT[["change_over_training"]] <- a_change$audit
