# F03 composites: the four within-group responses as the lead figure, the five
# between-responder contrasts as the second. Assembly only - the ring grids are built by
# the panel scripts. No statistics.
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

composite_responses <- annotate_composite(
  F03_PANELS$within_group,
  "F03 · How each arm responds: pathway enrichment on the volcano",
  paste(
    "Training (A-B) and acute (C-D) responses; HR (blue, left column) beside LR (red, right column).",
    "Both arms respond strongly; whether they respond differently is the between-responder test in the second figure."
  )
)

composite_contrasts <- annotate_composite(
  F03_PANELS$responders,
  "F03 · The between-responder test: HR vs LR states and response interactions",
  paste(
    "Between-responder states (A-C) and differential-response interactions (D-E).",
    "Arcs carry a genuine fgsea BH q; volcano points are pi-gated. No protein survives FDR in any of these contrasts."
  )
)
