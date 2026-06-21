# F01 Full Phenotype Atlas — Composite Figure
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

RPT <- "04_Figures/F01/b_reports"
pdf_device <- get_pdf_device()

read_panel_png <- function(filename) {
  path <- file.path(RPT, filename)
  if (!file.exists(path)) stop("Missing panel: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

wrap_panel <- function(grob) {
  ggplot() + annotation_custom(grob) +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
}

pA <- wrap_panel(read_panel_png("panel_A_volume_load_MAIN.png"))
pG <- wrap_panel(read_panel_png("panel_G_hypertrophy_composite_MAIN.png"))
pB <- wrap_panel(read_panel_png("panel_B_1rm_leg_MAIN.png"))
pH <- wrap_panel(read_panel_png("panel_H_1rm_ext_MAIN.png"))
pC <- wrap_panel(read_panel_png("panel_C_fcsa_mixed_MAIN.png"))
pI <- wrap_panel(read_panel_png("panel_I_mcsa_MAIN.png"))
pD <- wrap_panel(read_panel_png("panel_D_fcsa_type1_MAIN.png"))
pE <- wrap_panel(read_panel_png("panel_E_fcsa_type2_MAIN.png"))
pF <- wrap_panel(read_panel_png("panel_F_fiber_counts_MAIN.png"))

design <- "
AB
CD
EF
GH
II
"

composite <- pA + pG + pB + pH + pC + pI + pD + pE + pF +
  plot_layout(design = design, heights = c(0.65, 1, 1, 1, 1.5)) +
  plot_annotation(
    title = "Phenotype Atlas",
    subtitle = "Training load, hypertrophy, strength, muscle size, and fiber-resolved phenotypes",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

COMP_W <- 360
COMP_H <- 520

ggsave(file.path(RPT, "MAIN_F01_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "MAIN_F01_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("F01 composite saved")
