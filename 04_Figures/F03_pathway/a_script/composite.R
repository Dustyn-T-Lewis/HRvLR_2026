# F03 composite: the five responder contrasts as one figure. Assembly only - the ring
# grid is already built by panel_a_responders. No statistics.
pacman::p_load(ggplot2, patchwork)

composite <- F03_PANELS$responders +
  plot_annotation(
    title = "F03 · Where high and low responders differ: pathway enrichment on the volcano",
    subtitle = paste(
      "Between-responder contrasts (A-C) and the differential-response interactions (D-E).",
      "Arcs are fgsea BH q; points are pi-gated. Nothing survives FDR."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )
