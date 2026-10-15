library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F4")
stopifnot(names(panels) == c(
  "a_trait_changes", "b_links_vs_chance", "c_trait_pathways", "d_profile_clusters"
))
plots <- purrr::map2(lapply(panels, `[[`, "panel"), letters[seq_along(panels)], tag_panel)
# free() keeps panels a, b and d from aligning their axes to the long set names of panel c.
design <- c(area(1, 1, 1, 75), area(1, 76, 1, 180), area(2, 1, 2, 180), area(3, 1, 3, 180))
figure <- wrap_plots(
  free(plots$a_trait_changes), free(plots$b_links_vs_chance), plots$c_trait_pathways,
  free(plots$d_profile_clusters),
  design = design
) +
  plot_layout(heights = c(70, 90, 75))
ggsave(
  here("04_Figures", "F4", "b_reports", "F4_composite.pdf"), figure,
  width = 180, height = 240, units = "mm", device = cairo_pdf
)

write_figure_workbook(panels, "F4")
