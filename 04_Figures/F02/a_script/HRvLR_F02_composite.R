# F02 composite: the six main panels assembled live with patchwork so the plot panels
# align on a common grid (equal panel borders). Titles and subtitles are per-panel;
# tags come from patchwork; all legends collect to the figure bottom so no panel loses
# height to its own legend. Sources the panel scripts to get the ggplot objects pA..pF.
# Row 1: A state, B trajectory, C divergence; Row 2: D effect size, E DEPs, F pathways.
pacman::p_load(here, ggplot2, patchwork, grid)

source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))
for (p in c("A", "B", "C", "D", "E", "F")) {
  source(here("04_Figures", "F02", "a_script", "panels", sprintf("panel_%s.R", p)))
}

titled <- function(p, ttl, sub) {
  p + labs(title = ttl, subtitle = sub) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(face = "italic", size = 6, colour = "grey30"),
      plot.margin = margin(4, 4, 3, 3)
    )
}

grid <- titled(pA, "Global Proteome State", "PCA — no detectable group separation (PERMANOVA p = 0.62)") +
  titled(pB, "Response Magnitude ‖Δ‖", "RRPP — no detectable magnitude difference") +
  titled(pC, "Protein-level divergence leads", "pi-gated candidates; none survive FDR (exploratory)") +
  titled(pE, "DEPs per Contrast", "% of proteome (light = nominal p, dark = pi); dotted = chance") +
  titled(pD, "Effect-Size Distribution", "logFC density per contrast; number = median |log2FC|") +
  titled(pF, "Pathway Enrichment (GO)", "non-redundant sets, BH < 0.05") +
  plot_layout(ncol = 3, widths = c(1.3, 1, 1.35), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.margin = margin(2, 2, 2, 2))
  ) &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 11),
    plot.tag.position = c(0, 1)
  )

ggsave(file.path(RPT_DIR, "F02_composite.png"), grid,
  width = 310, height = 175, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT_DIR, "F02_composite.pdf"), grid,
  width = 235, height = 175, units = "mm", device = PDF_DEVICE, bg = "white"
)
message("F02 composite saved to ", RPT_DIR)
