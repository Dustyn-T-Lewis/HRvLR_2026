# S6: elastic net on traits, everyone pooled. Row 1: every fit. Row 2: what carries the lead, and
# whether it holds in both groups. Row 3: the lead's held-out predictions.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S6")
stopifnot(names(panels) == c("a_fit_overview", "b_selected", "c_predictions", "d_hr_lr"))
plots <- lapply(panels, \(env) free(env$panel))

# Widths in 5 mm columns of a 250 mm page: a across; b and d 125 mm each; c across the bottom left.
layout <- c(area(1, 1, 1, 50), area(2, 1, 2, 25), area(3, 1, 3, 25), area(2, 26, 3, 50))
figure <- wrap_plots(plots, design = layout, heights = c(80, 110, 70)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, "S6", 260)

write_figure_workbook(panels, "S6")
