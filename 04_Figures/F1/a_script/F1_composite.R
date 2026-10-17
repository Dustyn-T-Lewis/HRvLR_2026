# F1: design and the proteome in row 1, trait change, protein counts and fry transfer in row 2.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F1")
stopifnot(
  names(panels) == c("a_design", "b_pca", "c_trait_effects", "d_counts", "e_fry_transfer")
)
# Panels are not aligned with each other, so a wrapped panel's title stays beside its tag.
plots <- lapply(panels, \(env) free(env$panel))

# Widths in 5 mm columns: a 75 mm beside b 105 mm, then c 50, d 50 and e 80 mm.
layout <- c(
  area(1, 1, 1, 15), area(1, 16, 1, 36),
  area(2, 1, 2, 10), area(2, 11, 2, 20), area(2, 21, 2, 36)
)
figure <- wrap_plots(plots, design = layout, heights = c(78, 100)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
ggplot2::ggsave(
  here("04_Figures", "F1", "b_reports", "F1_composite.pdf"), figure,
  width = 180, height = 185, units = "mm", device = grDevices::cairo_pdf
)

write_figure_workbook(panels, "F1")
