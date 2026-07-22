# Panel B: every subject's composite score, ranked. Carries both the split (the
# HR/LR boundary) and its magnitude as one continuum, so the defining score needs
# only one panel. Circular by construction (the score defines the groups).
if (!exists("meta")) source(here::here("04_Figures", "F01_phenotype", "a_script", "setup.R"))
pacman::p_load(ggplot2, forcats)

build_continuum <- function(meta, tag = "B") {
  comp <- f01_composite_scores(meta) |>
    arrange(value) |>
    mutate(subject = fct_inorder(subject))
  boundary <- (min(comp$value[comp$Group == "HR"]) + max(comp$value[comp$Group == "LR"])) / 2
  swc <- 0.2 * sd(comp$value)

  p <- ggplot(comp, aes(value, subject, color = Group)) +
    annotate("rect",
      xmin = boundary - swc, xmax = boundary + swc, ymin = -Inf, ymax = Inf,
      fill = "grey85", alpha = 0.5
    ) +
    geom_vline(xintercept = boundary, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = value, yend = subject), linewidth = 0.5) +
    geom_point(size = 2.4) +
    scale_color_manual(values = GROUP_COLORS, name = NULL) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = "The responder split, as a continuum", tag = tag,
      subtitle = "Per-subject composite %change; dashed = HR/LR split, grey = small-difference zone",
      x = "Composite hypertrophy (%)", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40"),
      axis.text.y = element_text(size = 6.5),
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.2),
      panel.grid.major.y = element_blank()
    )
  audit <- comp |>
    transmute(subject, Group, composite = value, boundary, swc)
  list(plot = p, audit = audit)
}

b_continuum <- build_continuum(meta, tag = "B")
save_panel(b_continuum$plot, file.path(F01_RPT, "panels", "panel_b_continuum"), 150, 120)
F01_PANELS[["continuum"]] <- b_continuum$plot
F01_AUDIT[["responder_axis"]] <- b_continuum$audit
