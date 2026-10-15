library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F2")
stopifnot(names(panels) == c("a_responders", "b_hr_lr_candidates", "c_trajectories"))

# a and b share one colour cap, so their legends are identical and drawn once.
stopifnot(identical(
  panels$a_responders$panel$scales$get_scales("colour")$limits,
  panels$b_hr_lr_candidates$panel$scales$get_scales("colour")$limits
))
top_row <- (panels$a_responders$panel | panels$b_hr_lr_candidates$panel) +
  plot_layout(guides = "collect")
figure <- (top_row / panels$c_trajectories$panel) +
  plot_layout(heights = c(3.2, 1)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
ggsave(
  here("04_Figures", "F2", "b_reports", "F2_composite.pdf"), figure,
  width = 180, height = 215, units = "mm", device = cairo_pdf
)

write_figure_workbook(panels, "F2")
