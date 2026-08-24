# Panel B: the composite every subject was ranked on, and where the cut fell.
# Circular by construction - the composite defines the groups - so the panel
# reports the construction rather than treating the separation as a result.
if (!exists("ranked")) {
  source(here::here("04_Figures", "F01_phenotype", "a_script", "setup.R"))
}
pacman::p_load(ggplot2)

build_continuum <- function(ranked, composite_modality, tag = "B") {
  cut_at <- mean(c(
    min(ranked$comp_hypertrophy[ranked$group_arm == "HR"]),
    max(ranked$comp_hypertrophy[ranked$group_arm == "LR"])
  ))

  p <- ggplot(ranked, aes(.data$comp_hypertrophy, .data$subject,
    colour = .data$group_arm
  )) +
    geom_vline(
      xintercept = cut_at, linetype = "dashed", colour = "grey50",
      linewidth = 0.4
    ) +
    geom_segment(aes(
      x = 0, xend = .data$comp_hypertrophy,
      yend = .data$subject
    ), linewidth = 0.5) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = GROUP_COLORS, name = NULL) +
    labs(
      title = "The split, and the score it was cut from", tag = tag,
      subtitle = sprintf(
        "Dashed = HR/LR cut. Composite is bimodal: bootstrap LRT p = %.3f",
        composite_modality$lrt_p[1]
      ),
      x = "Composite hypertrophy", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      axis.text.y = element_text(size = 6.5),
      legend.position = "inside",
      legend.position.inside = c(0.86, 0.18),
      panel.grid.major.y = element_blank()
    )

  list(
    plot = p,
    audit = ranked |>
      transmute(subject, group_arm, comp_hypertrophy, cut_at)
  )
}

b_continuum <- build_continuum(ranked, composite_modality)
save_panel(
  b_continuum$plot, file.path(F01_RPT, "panels", "panel_b_continuum"),
  135, 120
)
F01_PANELS[["continuum"]] <- b_continuum$plot
F01_AUDIT[["responder_axis"]] <- b_continuum$audit
