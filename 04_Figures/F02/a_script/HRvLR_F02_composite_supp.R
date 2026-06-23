# F02 supplementary composite — variability panels + pilot p < 0.05 Venns.
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

composite <- (panel("panel_g_upset") | panel("panel_n_venn_primary")) /
  (panel("panel_c_cv") | panel("panel_d_cv_transitions") | panel("panel_e_intra_variability")) /
  (panel("panel_k_venn_training_p") | panel("panel_l_venn_acute_p") | plot_spacer()) +
  plot_layout(heights = c(1, 1, 0.9)) +
  plot_annotation(
    title = "Global Proteome Overview — Supplementary",
    subtitle = "Contrast overlap (UpSet, 5-set Euler), variability, and nominal-p HR-vs-LR overlap",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT, "F02_composite_supp.png"), composite,
  width = 420, height = 380, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT, "F02_composite_supp.pdf"), composite,
  width = 420, height = 380, units = "mm", device = PDF_DEVICE, bg = "white"
)
cat("F02 supplementary composite saved to", RPT, "\n")
