# S5: module detail. Row 1: what each module holds, and the soft-power QC beside it. Row 2: STRING
# edges within each module, and the module with the smallest eigengene p person by person.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S5")
stopifnot(names(panels) == c("a_characterisation", "b_string", "c_module_qc", "d_trajectories"))
plots <- purrr::map2(lapply(panels, `[[`, "panel"), letters[seq_along(panels)], tag_panel)

# Widths on a 250 mm page: a 155 mm beside c 95 mm; then b 140 mm beside d 110 mm. Each panel is
# wrapped so its axes do not line up with the others'.
row <- \(left, right, widths) {
  wrap_elements(full = wrap_plots(
    wrap_elements(full = left), wrap_elements(full = right),
    widths = widths
  ))
}
figure <- wrap_plots(
  row(plots$a_characterisation, plots$c_module_qc, c(155, 95)),
  row(plots$b_string, plots$d_trajectories, c(140, 110)),
  ncol = 1, heights = c(115, 95)
)
save_figure(figure, "S5", 210)

write_figure_workbook(panels, "S5")
