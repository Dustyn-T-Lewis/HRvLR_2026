# F1: design and the proteome in row 1; trait change and protein and gene-set counts in row 2; fry
# transfer in row 3.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F1")
stopifnot(
  names(panels) == c(
    "a_design", "b_pca", "c_trait_effects", "d_counts", "e_set_counts", "f_fry_transfer"
  )
)
# Panels are not aligned with each other, so a wrapped panel's title stays beside its tag.
plots <- lapply(panels, \(env) free(env$panel))

# Widths in 5 mm columns of a 250 mm page: a 95 mm beside b 155 mm; then c 80, d 90 and e 80 mm,
# so proteins (d) and gene sets (e) per contrast sit side by side; f across the page.
layout <- c(
  area(1, 1, 1, 19), area(1, 20, 1, 50),
  area(2, 1, 2, 16), area(2, 17, 2, 34), area(2, 35, 2, 50), area(3, 1, 3, 50)
)
figure <- wrap_plots(plots, design = layout, heights = c(66, 95, 70)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, "F1", 235)

write_figure_workbook(panels, "F1")
