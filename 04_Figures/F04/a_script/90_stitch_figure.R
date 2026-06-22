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

# Nine panels grouped into four family rows: HR, LR, HRvLR, Interaction.
pA <- read_panel_png(main_dir, "MAIN_panel_A_Training_HR.png")
pB <- read_panel_png(main_dir, "MAIN_panel_B_Acute_HR.png")
pC <- read_panel_png(main_dir, "MAIN_panel_C_Training_LR.png")
pD <- read_panel_png(main_dir, "MAIN_panel_D_Acute_LR.png")
pE <- read_panel_png(main_dir, "MAIN_panel_E_Baseline_HRvLR.png")
pF <- read_panel_png(main_dir, "MAIN_panel_F_Trained_HRvLR.png")
pG <- read_panel_png(main_dir, "MAIN_panel_G_Acute_HRvLR.png")
pH <- read_panel_png(main_dir, "MAIN_panel_H_Training_Interaction.png")
pI <- read_panel_png(main_dir, "MAIN_panel_I_Acute_Interaction.png")

wrap_panel <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

main_composite <- (
  wrap_panel(pA) | wrap_panel(pB) | plot_spacer()
) / (
  wrap_panel(pC) | wrap_panel(pD) | plot_spacer()
) / (
  wrap_panel(pE) | wrap_panel(pF) | wrap_panel(pG)
) / (
  wrap_panel(pH) | wrap_panel(pI) | plot_spacer()
) +
  plot_annotation(
    title = "Protein- and Pathway-Level Signatures",
    subtitle = "Volcano rings by family: HR and LR within-responder, HR-vs-LR state, and differential response",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, face = "bold.italic", color = "grey30")
    )
  )

dir.create(file.path(RPT_MAIN, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_MAIN, "png"), recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(RPT_MAIN, "pdf", "MAIN_F04_composite.pdf"), main_composite,
  width = 540, height = 740, units = "mm",
  device = pdf_device, limitsize = FALSE
)
ggsave(file.path(RPT_MAIN, "png", "MAIN_F04_composite.png"), main_composite,
  width = 540, height = 740, units = "mm",
  dpi = 300, limitsize = FALSE
)

cat("F04 composite saved\n")
