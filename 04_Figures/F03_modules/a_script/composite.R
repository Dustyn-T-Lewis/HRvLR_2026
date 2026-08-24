# F03 composite: the atlas and the calibration on the left, the full grid on
# the right. Assembly only, no statistics.
pacman::p_load(patchwork, ggplot2)

left_col <- (F03_PANELS$atlas / F03_PANELS$calibration) +
  plot_layout(heights = c(1.4, 1))

composite <- (left_col | F03_PANELS$heatmap) +
  plot_layout(widths = c(1, 1.35)) +
  plot_annotation(
    title = "F03 · Co-expression modules against adaptation",
    subtitle = paste(
      "Twelve modules named by GO, tested against ten phenotypes in six",
      "windows, and read against a shuffled-phenotype null."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )
  )
