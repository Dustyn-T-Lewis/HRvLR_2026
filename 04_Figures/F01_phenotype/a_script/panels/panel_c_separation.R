# Panel C: what the label separates, on the same scale as panel A so the two
# read together. Point shape marks whether the outcome is an ingredient of the
# composite the label was cut from, because an ingredient cannot corroborate it.
if (!exists("label_separation")) {
  source(here::here("04_Figures", "F01_phenotype", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats)

build_separation <- function(label_separation, composite_structure, tag = "C") {
  d <- label_separation |>
    filter(.data$trait != "comp_hypertrophy") |>
    left_join(composite_structure, by = "trait") |>
    mutate(
      label = TRAIT_LABELS[.data$trait],
      internal = .data$r2_alone >= 0.20,
      separated = .data$ci_lo_d > 0 | .data$ci_hi_d < 0,
      label = factor(.data$label, levels = TRAIT_ORDER)
    )

  p <- ggplot(d, aes(.data$difference_d, .data$label,
    colour = .data$separated,
    shape = .data$internal
  )) +
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
    scale_shape_manual(
      values = c(`TRUE` = 17, `FALSE` = 16),
      labels = c(`TRUE` = "in composite", `FALSE` = "independent"),
      name = NULL
    ) +
    labs(
      title = "What the HR/LR label separates", tag = tag,
      subtitle = "HR minus LR, in SD units, with 95% CI",
      x = "Group difference (SD units)", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      legend.position = "inside",
      legend.position.inside = c(0.17, 0.17),
      legend.background = element_rect(
        fill = alpha("white", 0.7),
        colour = NA
      ),
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = d |> select(
    trait, hr, lr, difference, ci_lo, ci_hi,
    difference_d, ci_lo_d, ci_hi_d, p, r2_alone, internal, separated
  ))
}

c_separation <- build_separation(label_separation, composite_structure)
save_panel(
  c_separation$plot,
  file.path(F01_RPT, "panels", "panel_c_separation"), 135, 85
)
F01_PANELS[["separation"]] <- c_separation$plot
F01_AUDIT[["label_separation"]] <- c_separation$audit
