# F02 composite: the whole sweep across the top, the one survivor and its
# robustness beneath. Assembly only, no statistics.
pacman::p_load(patchwork, ggplot2)

bottom <- F02_PANELS$hit | F02_PANELS$robust

composite <- F02_PANELS$landscape / bottom +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = "F02 · Does a proteome change track an adaptation?",
    subtitle = paste(
      "Per-subject change regressed on per-subject adaptation.",
      "No groups, no cut points, no baseline."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )
  )
