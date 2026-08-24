# F02 composite: the null comparison on top because it carries the result, the
# false-positive rate and the agreement below it. Assembly only, no statistics.
pacman::p_load(patchwork, ggplot2)

composite <- F02_PANELS$null /
  (F02_PANELS$falsepos | F02_PANELS$agreement) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "F02 · Does the proteome group these subjects on its own?",
    subtitle = paste(
      "Baseline proteome clustered blind to every label, calibrated against",
      "data with no cluster structure."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )
  )
