# F2: responders and HR - LR candidates as dot-matrices sharing one key, then the candidates with
# the smallest p person by person.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F2")
stopifnot(names(panels) == c("a_responders", "b_hr_lr_candidates", "c_trajectories"))

# Widths on a 250 mm page: a and b 80 mm each, their shared key 18 mm, c 72 mm. Both matrices take
# protein_colour_limit() of the same table, so their keys are identical and collected once.
# c keeps its own key at its foot.
figure <- wrap_plots(
  panels$a_responders$panel, panels$b_hr_lr_candidates$panel, guide_area(),
  panels$c_trajectories$panel + plot_layout(guides = "keep"),
  widths = c(80, 80, 18, 72), guides = "collect"
) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, "F2", 140)

write_figure_workbook(panels, "F2")
