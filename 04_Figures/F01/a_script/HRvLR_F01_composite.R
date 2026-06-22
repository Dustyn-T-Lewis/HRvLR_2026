# F01 composite — stitch the rendered phenotype panels into one figure.
# Run after HRvLR_F01_run.R has written the panels to b_reports.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

RPT <- "04_Figures/F01/b_reports"

panel <- function(name) {
  path <- file.path(RPT, paste0(name, ".png"))
  if (!file.exists(path)) stop("Missing panel: ", path, " (run HRvLR_F01_run.R first)")
  ggplot() +
    annotation_custom(rasterGrob(readPNG(path), interpolate = TRUE)) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

panels <- c(
  "panel_a_volume_load", "panel_b_1rm_leg", "panel_c_fcsa_mixed",
  "panel_d_fcsa_type1", "panel_e_fcsa_type2", "panel_f_fiber_counts",
  "panel_g_hypertrophy_composite", "panel_h_1rm_ext", "panel_i_mcsa"
)

composite <- wrap_plots(lapply(panels, panel), ncol = 3) +
  plot_annotation(
    title = "Phenotype Atlas (HR vs LR)",
    subtitle = "Training volume, strength, and fiber/muscle cross-sectional area",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT, "F01_composite.png"), composite,
  width = 330, height = 330, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT, "F01_composite.pdf"), composite,
  width = 330, height = 330, units = "mm", device = PDF_DEVICE, bg = "white"
)
cat("F01 composite saved to", RPT, "\n")
