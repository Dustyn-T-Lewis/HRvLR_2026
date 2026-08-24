# F01 composite: what adapted on the left over the subject map, the
# redundancy structure on the right. Assembly only, no statistics.
pacman::p_load(patchwork, ggplot2)

left_col <- (F01_PANELS$change / F01_PANELS$space) +
  plot_layout(heights = c(1, 1.25))

composite <- (left_col | F01_PANELS$structure) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "F01 · The adaptations this proteome is mapped onto",
    subtitle = paste(
      "Ten phenotypes treated as continuous outcomes: what moved, how the",
      "subjects sit, and how much the measures overlap."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )
  )
