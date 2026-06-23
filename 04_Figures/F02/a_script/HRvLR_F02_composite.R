# F02 composite — stitch the rendered proteome-overview panels into one figure.
# Run after HRvLR_F02_run.R has written the panels to b_reports.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

RPT <- "04_Figures/F02/b_reports"

panel <- function(name) {
  path <- file.path(RPT, paste0(name, ".png"))
  if (!file.exists(path)) stop("Missing panel: ", path, " (run HRvLR_F02_run.R first)")
  ggplot() +
    annotation_custom(rasterGrob(readPNG(path), interpolate = TRUE)) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# Row 2: DEP bars + histograms + pathways, equal height, histograms a touch narrower
row_dab <- panel("panel_f_dep_counts") + panel("panel_b_effect_size") +
  panel("panel_h_pathways") + plot_layout(widths = c(1, 0.82, 1))
# Row 3: the two Venns, smaller (flanked by spacers)
row_venn <- plot_spacer() + panel("panel_i_venn_training") +
  panel("panel_j_venn_acute") + plot_spacer() + plot_layout(widths = c(0.35, 1, 1, 0.35))

composite <- panel("panel_a_pca") / row_dab / row_venn / panel("panel_m_eta2") +
  plot_layout(heights = c(1.05, 1, 0.8, 0.55)) +
  plot_annotation(
    title = "Global Proteome Overview",
    subtitle = "Structure, differential abundance, contrast overlap, and design-structured variance",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT, "F02_composite.png"), composite,
  width = 400, height = 440, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT, "F02_composite.pdf"), composite,
  width = 400, height = 440, units = "mm", device = PDF_DEVICE, bg = "white"
)
cat("F02 composite saved to", RPT, "\n")
