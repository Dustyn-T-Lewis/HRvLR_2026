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

# Widths in 5 mm columns of a 250 mm page: a 95 mm beside b 155 mm, then c 70, d 75 and e 105 mm.
layout <- c(
  area(1, 1, 1, 19), area(1, 20, 1, 50),
  area(2, 1, 2, 14), area(2, 15, 2, 29), area(2, 30, 2, 50)
)
figure <- wrap_plots(plots, design = layout, heights = c(66, 95)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, "F1", 165)

write_figure_workbook(panels, "F1")
