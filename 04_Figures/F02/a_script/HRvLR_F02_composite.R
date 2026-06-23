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

composite <- panel("panel_a_pca") /
  (panel("panel_f_dep_counts") | panel("panel_b_effect_size") | panel("panel_h_pathways")) /
  (panel("panel_g_upset") | panel("panel_i_venn_training") | panel("panel_j_venn_acute")) /
  (panel("panel_c_cv") | panel("panel_d_cv_transitions") | panel("panel_e_intra_variability")) +
  plot_layout(heights = c(0.8, 1, 0.9, 0.9)) +
  plot_annotation(
    title = "Global Proteome Overview",
    subtitle = "Structure, differential abundance, contrast overlap, and variability",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT, "F02_composite.png"), composite,
  width = 420, height = 480, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT, "F02_composite.pdf"), composite,
  width = 420, height = 480, units = "mm", device = PDF_DEVICE, bg = "white"
)
cat("F02 composite saved to", RPT, "\n")
