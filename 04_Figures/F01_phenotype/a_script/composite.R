# F01 composite: the matched-work control and the responder continuum stack on
# the left, the divergence forest carries the right. Assembly only, no statistics.
pacman::p_load(patchwork, ggplot2)

left_col <- (F01_PANELS$volume / F01_PANELS$continuum) +
  plot_layout(heights = c(1, 1.9))

composite <- (left_col | F01_PANELS$forest) +
  plot_layout(widths = c(1, 1.1)) +
  plot_annotation(
    title = "F01 · Phenotype atlas: matched work, divergent growth",
    subtitle = "Matched training work; groups classified by growth magnitude; divergence by domain.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )
