# F02 composite: the six main panels assembled with patchwork so the plot panels align
# on a common grid. Titles and subtitles are per-panel; tags come from patchwork; all
# legends collect to the figure bottom so no panel loses height to its own legend.
# Row 1: A state, B trajectory, C divergence; Row 2: D effect size, E DEPs, F pathways.
# Assembly only - the panels (pA..pF) are already in scope from 01_run_proteome.R. No statistics.
pacman::p_load(ggplot2, patchwork, grid)

titled <- function(p, ttl, sub) {
  p + labs(title = ttl, subtitle = sub) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(face = "italic", size = 6, colour = "grey30"),
      plot.margin = margin(4, 4, 3, 3)
    )
}

composite <- titled(pA, "Global Proteome State", sprintf(
  "PCA — no detectable group separation (PERMANOVA p = %.2f)", perm_group[["Pr(>F)"]][1]
)) +
  titled(pB, "Response Magnitude ||delta||", "RRPP — no detectable magnitude difference") +
  titled(pC, "Protein-level divergence leads", "pi-gated candidates; none survive FDR (exploratory)") +
  titled(pD, "DEPs per Contrast", "% of proteome (light = nominal p, dark = pi); dotted = chance") +
  titled(pE, "Effect-Size Distribution", "logFC density per contrast; number = median |log2FC|") +
  titled(pF, "Pathway Enrichment (GO)", "non-redundant sets, BH < 0.05") +
  plot_layout(ncol = 3, widths = c(1.3, 1, 1.35)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.margin = margin(2, 2, 2, 2))
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 11),
    plot.tag.position = c(0, 1)
  )
