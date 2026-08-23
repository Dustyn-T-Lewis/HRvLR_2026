# F03 composite, continuous tree: one figure, the two pooled contrasts.
# Mirrors categorical/F03_pathway/a_script/composite.R's
# annotate_composite() pattern -- one figure instead of two, since there's
# no between-arm contrast set to give a second figure to. Assembly only,
# no statistics.
pacman::p_load(ggplot2, patchwork)

annotate_composite <- function(grid, title, subtitle) {
  grid +
    plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey30")
      )
    )
}

composite_pooled <- annotate_composite(
  F03_PANELS$pooled,
  "F03 (continuous) - Pathway enrichment on the volcano, no arm split",
  paste(
    "Training (A) and acute (B) responses, pooled across all 16 subjects.",
    "Arcs carry a genuine fgsea BH q; volcano points are pi-gated."
  )
)
