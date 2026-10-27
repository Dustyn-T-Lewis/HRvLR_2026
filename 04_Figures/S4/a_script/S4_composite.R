# S4: slope-difference counts across the page, then the strongest screen and its top feature.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S4")
stopifnot(names(panels) == c("a_slope_counts", "b_within_rho", "c_top_slope"))

# Widths in 5 mm columns of a 250 mm page: a full width, then b and c 125 mm each.
layout <- c(area(1, 1, 1, 50), area(2, 1, 2, 25), area(2, 26, 2, 50))
figure <- wrap_plots(lapply(panels, \(env) free(env$panel)), design = layout, heights = c(65, 95)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, "S4", 160)

write_figure_workbook(panels, "S4")
