# F06 HRvLR — Composite: WGCNA Module Analysis
# Layout: A (top) / (B | C) (bottom)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

RPT <- "04_Figures/F06/b_reports"
pdf_device <- get_pdf_device()

read_panel_png <- function(filename) {
  path <- file.path(RPT, filename)
  if (!file.exists(path)) stop("Missing panel: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA <- read_panel_png("panel_A_module_trait_MAIN.png")
pB <- read_panel_png("panel_B_triptych_MAIN.png")
pC <- read_panel_png("panel_C_hub_network_MAIN.png")

wrap_panel <- function(grob) {
  ggplot() + annotation_custom(grob) +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
}

composite <- wrap_panel(pA) /
             (wrap_panel(pB) | wrap_panel(pC)) +
  plot_layout(heights = c(0.35, 0.65))

COMP_W <- 400; COMP_H <- 600

ggsave(file.path(RPT, "F06_wgcna_MAIN.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "F06_wgcna_MAIN.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("F06 HRvLR WGCNA composite saved\n")
