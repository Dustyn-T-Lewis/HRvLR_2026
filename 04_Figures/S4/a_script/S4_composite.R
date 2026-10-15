library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S4")
stopifnot(names(panels) == c("a_slope_counts", "b_within_rho", "c_top_slope"))

figure <- panels$a_slope_counts$panel /
  (panels$b_within_rho$panel | panels$c_top_slope$panel) +
  plot_layout(heights = c(65, 90)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
ggsave(
  here("04_Figures", "S4", "b_reports", "S4_composite.pdf"), figure,
  width = 180, height = 165, units = "mm", device = cairo_pdf
)

write_figure_workbook(panels, "S4")
