# F01 composite, continuous tree: the composite continuum on the left, the
# change-magnitude summary on the right. Two panels, not three -- the
# categorical tree's matched-training-volume panel has no group-free
# analog (there is nothing to match without two arms), so it is dropped
# rather than forced.
pacman::p_load(patchwork, ggplot2)

composite <- (F01_PANELS$continuum | F01_PANELS$magnitude) +
  plot_layout(widths = c(1, 1.1)) +
  plot_annotation(
    title = "F01 (continuous) · Phenotype: adaptation as a continuum",
    subtitle = paste0(
      "No group split anywhere: every subject ranked by composite score, ",
      "every outcome's change read as a distribution."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )
