library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S5")
stopifnot(names(panels) == c("a_characterisation", "b_string", "c_module_qc"))
plots <- purrr::map2(lapply(panels, `[[`, "panel"), letters[seq_along(panels)], tag_panel)
# Each panel is wrapped so its axes do not line up with the other panels' axes.
figure <- wrap_plots(lapply(plots, \(p) wrap_elements(full = p)),
  ncol = 1, heights = c(100, 85, 42)
)
ggsave(
  here("04_Figures", "S5", "b_reports", "S5_composite.pdf"), figure,
  width = 180, height = 238, units = "mm", device = cairo_pdf
)

write_figure_workbook(panels, "S5")
