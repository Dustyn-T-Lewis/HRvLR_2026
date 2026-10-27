# S5: module detail. Row 1: what each module holds, and the soft-power QC beside it. Row 2: STRING
# edges within each module.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S5")
stopifnot(names(panels) == c("a_characterisation", "b_string", "c_module_qc"))
plots <- purrr::map2(lapply(panels, `[[`, "panel"), letters[seq_along(panels)], tag_panel)

# Widths on a 250 mm page: a 155 mm beside c 95 mm, then b full width. Each panel is wrapped so
# its axes do not line up with the others'.
figure <- wrap_plots(
  wrap_elements(full = wrap_plots(
    wrap_elements(full = plots$a_characterisation), wrap_elements(full = plots$c_module_qc),
    widths = c(155, 95)
  )),
  wrap_elements(full = plots$b_string),
  ncol = 1, heights = c(115, 80)
)
save_figure(figure, "S5", 195)

write_figure_workbook(panels, "S5")
