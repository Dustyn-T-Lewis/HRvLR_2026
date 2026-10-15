library(here)
library(dplyr)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F1")
stopifnot(names(panels) == c("a_design", "b_pca", "c_counts", "d_fry_transfer"))
# Panels are not aligned with each other, so a wrapped panel's title stays beside its tag.
plots <- lapply(panels, \(env) free(env$panel))

# Widths in 5 mm columns: a 95 mm beside b 85 mm, then c 85 mm beside d 95 mm.
layout <- c(area(1, 1, 1, 19), area(1, 20, 1, 36), area(2, 1, 2, 17), area(2, 18, 2, 36))
figure <- wrap_plots(plots, design = layout, heights = c(88, 84)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
ggplot2::ggsave(
  here("04_Figures", "F1", "b_reports", "F1_composite.pdf"), figure,
  width = 180, height = 172, units = "mm", device = cairo_pdf
)

write_figure_workbook(panels, "F1")
