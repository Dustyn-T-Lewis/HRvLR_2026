# F5: modules. Row 1: every module across the contrasts and the trait links. Row 2: how modules
# change together after the bout within each group, and preservation between the groups.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F5")
stopifnot(names(panels) == c("a_heatmap", "b_coupling", "c_preservation"))
plots <- purrr::map2(lapply(panels, `[[`, "panel"), letters[seq_along(panels)], tag_panel)

# Widths on a 250 mm page: a across the page; b 140 mm beside c 110 mm.
figure <- wrap_plots(
  wrap_elements(full = plots$a_heatmap),
  wrap_elements(
    full = wrap_plots(plots$b_coupling, plots$c_preservation, widths = c(140, 110))
  ),
  ncol = 1, heights = c(95, 90)
)
save_figure(figure, "F5", 185)

write_figure_workbook(panels, "F5")
