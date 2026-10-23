# F3: the pathway matrix, with the gene-set counts per contrast beside it.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F3")
stopifnot(names(panels) == c("a_pathways", "b_set_counts"))

# Widths in 5 mm columns of a 250 mm page: a 170 mm, b 80 mm in the top 110 mm, open below it.
height <- panels$a_pathways$panel_height
stopifnot(height <= 240)
layout <- c(area(1, 1, 2, 34), area(1, 35, 1, 50))
figure <- wrap_plots(lapply(panels, \(env) free(env$panel)),
  design = layout, heights = c(110, height - 110)
) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, "F3", height)

write_figure_workbook(panels, "F3")
