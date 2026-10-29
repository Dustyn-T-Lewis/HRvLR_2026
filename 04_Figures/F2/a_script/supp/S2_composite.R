# S2: can single features tell HR from LR? Row 1: the top features by AUC per task. Row 2: the
# exploratory HR - LR protein candidates and the four with the smallest p, person by person.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S2")
stopifnot(names(panels) == c("a_auc_top", "b_hr_lr_candidates", "c_trajectories"))
plots <- lapply(panels, \(env) free(env$panel))

# Widths in 5 mm columns of a 250 mm page: a full width; then b 160 mm (with its key) and c 90 mm.
layout <- c(area(1, 1, 1, 50), area(2, 1, 2, 32), area(2, 33, 2, 50))
figure <- wrap_plots(plots, design = layout, heights = c(105, 135)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, "S2", 240)

write_figure_workbook(panels, "S2")
