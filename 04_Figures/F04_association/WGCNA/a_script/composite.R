# F04 composite: the atlas and what the modules are, then the two null read-outs.
# Assembly only, no statistics.
pacman::p_load(ggplot2, patchwork)

composite <- (F04_PANELS$atlas / F04_PANELS$pathways / (F04_PANELS$fry | F04_PANELS$phenotype)) +
  plot_layout(heights = c(1, 1.5, 1.2)) +
  plot_annotation(
    title = "F04 · Co-expression modules: coherent biology, no responder signal",
    subtitle = paste0(
      "Modules are defined on within-subject co-regulation and scored on raw abundance.\n",
      "Biologically interpretable (B), but none moves with responder status (C) ",
      "and none predicts the size of the adaptation (D)."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )
