# F04 HRvLR — Main and supplementary volcano-ring composites
here::i_am(".here")
source(here::here("04_Figures/F04/a_script/style.R"))

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

RPT_MAIN <- here::here("04_Figures/F04/b_reports/main")
RPT_SUPP <- here::here("04_Figures/F04/b_reports/supp")
pdf_device <- get_pdf_device()

read_panel_png <- function(dir_path, filename) {
  path <- file.path(dir_path, filename)
  if (!file.exists(path)) stop("Missing panel: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

main_dir <- file.path(RPT_MAIN, "png")
supp_dir <- file.path(RPT_SUPP, "png")

pA <- read_panel_png(main_dir, "MAIN_panel_A_Baseline_HRvLR.png")
pB <- read_panel_png(main_dir, "MAIN_panel_B_Trained_HRvLR.png")
pC <- read_panel_png(main_dir, "MAIN_panel_C_Acute_HRvLR.png")
pD <- read_panel_png(main_dir, "MAIN_panel_D_Training_HR.png")
pE <- read_panel_png(main_dir, "MAIN_panel_E_Training_LR.png")
pF <- read_panel_png(main_dir, "MAIN_panel_F_Acute_HR.png")
pG <- read_panel_png(main_dir, "MAIN_panel_G_Acute_LR.png")

sA <- read_panel_png(supp_dir, "SUPP_panel_A_Training_Interaction.png")
sB <- read_panel_png(supp_dir, "SUPP_panel_B_Acute_Interaction.png")

wrap_panel <- function(grob) {
  ggplot() + annotation_custom(grob) +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
}

main_composite <- (
  wrap_panel(pA) | wrap_panel(pB) | wrap_panel(pC)
) / (
  wrap_panel(pD) | wrap_panel(pE) | plot_spacer()
) / (
  wrap_panel(pF) | wrap_panel(pG) | plot_spacer()
) +
  plot_annotation(
    title = "Protein- and Pathway-Level Signatures",
    subtitle = "State contrasts plus within-group training and acute responses",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, face = "bold.italic", color = "grey30")
    )
  )

supp_composite <- (wrap_panel(sA) | wrap_panel(sB)) +
  plot_annotation(
    title = "Interaction Signatures",
    subtitle = "Difference-of-differences contrasts for training and acute response",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, face = "bold.italic", color = "grey30")
    )
  )

dir.create(file.path(RPT_MAIN, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_MAIN, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_SUPP, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_SUPP, "png"), recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(RPT_MAIN, "pdf", "MAIN_F04_composite.pdf"), main_composite,
       width = 540, height = 560, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_MAIN, "png", "MAIN_F04_composite.png"), main_composite,
       width = 540, height = 560, units = "mm",
       dpi = 300, limitsize = FALSE)

ggsave(file.path(RPT_SUPP, "pdf", "SUPP_F04_interactions.pdf"), supp_composite,
       width = 380, height = 210, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_SUPP, "png", "SUPP_F04_interactions.png"), supp_composite,
       width = 380, height = 210, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("F04 composites saved\n")
