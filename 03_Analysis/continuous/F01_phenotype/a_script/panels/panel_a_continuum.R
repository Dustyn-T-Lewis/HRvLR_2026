# Panel A: every subject's composite score, ranked. Mirrors
# categorical/F01_phenotype/a_script/panels/panel_b_continuum.R with the
# HR/LR boundary and the responder colour removed -- this tree never draws
# that line. What's left is the plain question: how much did people
# actually adapt, and how continuous is that spread.
if (!exists("meta")) {
  source(here::here(
    "03_Analysis", "continuous", "F01_phenotype", "a_script", "setup.R"
  ))
}
pacman::p_load(ggplot2, forcats)

build_continuum <- function(meta, tag = "A") {
  comp <- f01_composite_scores(meta) |>
    arrange(value) |>
    mutate(subject = fct_inorder(subject))

  p <- ggplot(comp, aes(value, subject)) +
    geom_segment(aes(x = 0, xend = value, yend = subject), linewidth = 0.5) +
    geom_point(size = 2.4, color = DOMAIN_COLORS[["muscle"]]) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = "Composite hypertrophy, no split", tag = tag,
      subtitle =
        "Every subject's composite %change, ranked; no group boundary drawn",
      x = "Composite hypertrophy (%)", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", color = "grey40"),
      axis.text.y = element_text(size = 6.5),
      panel.grid.major.y = element_blank()
    )
  list(plot = p, audit = transmute(comp, subject, composite = value))
}

a_continuum <- build_continuum(meta, tag = "A")
save_panel(
  a_continuum$plot, file.path(F01_RPT, "panels", "panel_a_continuum"), 150,
  120
)
F01_PANELS[["continuum"]] <- a_continuum$plot
F01_AUDIT[["composite_axis"]] <- a_continuum$audit
