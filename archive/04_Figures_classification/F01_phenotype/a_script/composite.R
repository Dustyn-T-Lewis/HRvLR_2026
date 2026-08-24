# F01 composite: the two effect-size panels stack on the left so they share an
# axis and read against each other; the continuum carries the right. Assembly
# only, no statistics.
pacman::p_load(patchwork, ggplot2)

left_col <- (F01_PANELS$change / F01_PANELS$separation) +
  plot_layout(heights = c(1, 1))

composite <- (left_col | F01_PANELS$continuum) +
  plot_layout(widths = c(1, 0.95)) +
  plot_annotation(
    title = "F01 · The responder phenotype and the label cut from it",
    subtitle = paste(
      "What changed over training, what the HR/LR label separates,",
      "and the composite the split was drawn on."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )
  )
