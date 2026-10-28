# F5: modules. Row 1: every module across the contrasts and the trait links. Row 2: the module with
# the smallest p person by person, and preservation between the groups.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F5")
stopifnot(names(panels) == c("a_heatmap", "b_trajectories", "c_preservation"))
plots <- purrr::map2(lapply(panels, `[[`, "panel"), letters[seq_along(panels)], tag_panel)

# Widths on a 250 mm page: a across the page; b 100 mm beside c 150 mm.
figure <- wrap_plots(
  wrap_elements(full = plots$a_heatmap),
  wrap_elements(
    full = wrap_plots(plots$b_trajectories, plots$c_preservation, widths = c(100, 150))
  ),
  ncol = 1, heights = c(95, 85)
)
save_figure(figure, "F5", 180)

write_figure_workbook(panels, "F5")
