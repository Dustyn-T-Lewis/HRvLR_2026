# F02 composite, continuous tree: three panels, assembled with patchwork.
# Mirrors categorical/F02_proteome/a_script/composite.R's per-panel-title
# pattern on a reduced panel set -- see the tree's own panel scripts for
# which categorical panels were dropped and why.
pacman::p_load(ggplot2, patchwork, grid)

titled <- function(p, ttl, sub) {
  p + labs(title = ttl, subtitle = sub) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(
        face = "italic", size = 6, colour = "grey30"
      ),
      plot.margin = margin(4, 4, 3, 3)
    )
}

composite <- titled(pA, "Global Proteome State", sprintf(
  "PCA -- hypertrophy PERMANOVA p = %.2f",
  perm_pheno[["Pr(>F)"]][1]
)) +
  titled(
    pB, "DEPs per Pooled Contrast",
    "% of proteome (light = nominal p, dark = pi); dotted = chance"
  ) +
  titled(
    pC, "Training vs acute concordance",
    "pi-gated both phases; none survive FDR (exploratory)"
  ) +
  plot_layout(ncol = 3, widths = c(1.3, 1, 1.2)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.margin = margin(2, 2, 2, 2))
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 11),
    plot.tag.position = c(0, 1)
  )
