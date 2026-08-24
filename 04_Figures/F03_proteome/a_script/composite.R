# F03 composite: the sweep on the left, its two calibration checks stacked on
# the right. Assembly only, no statistics.
pacman::p_load(patchwork, ggplot2)

right_col <- F03_PANELS$confirm / F03_PANELS$fgsea

composite <- (F03_PANELS$sweep | right_col) +
  plot_layout(widths = c(1, 1.05)) +
  plot_annotation(
    title = "F03 · What the sweep found, and whether it holds",
    subtitle = paste(
      "Two contrasts under six candidate labels, then the same tests run on",
      "randomly relabelled subjects."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )
  )
