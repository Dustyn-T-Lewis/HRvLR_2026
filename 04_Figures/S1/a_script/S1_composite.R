library(here)
library(dplyr)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S1")
stopifnot(names(panels) == c(
  "a_no_clusters", "b_variance_shares", "c_diffvar", "d_pc_covariates", "e_pca_levels"
))
# Panels are not aligned with each other, so a wrapped panel's title stays beside its tag.
plots <- lapply(panels, \(env) free(env$panel))

# Widths in 5 mm columns: a, c and e 110 mm, b and d 70 mm.
layout <- c(
  area(1, 1, 1, 22), area(1, 23, 1, 36), area(2, 1, 2, 22), area(2, 23, 2, 36),
  area(3, 1, 3, 22)
)
figure <- wrap_plots(plots, design = layout, heights = c(95, 72, 62)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
ggplot2::ggsave(
  here("04_Figures", "S1", "b_reports", "S1_composite.pdf"), figure,
  width = 180, height = 229, units = "mm", device = cairo_pdf
)

write_figure_workbook(panels, "S1")
