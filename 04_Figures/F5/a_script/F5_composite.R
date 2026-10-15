library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F5")
stopifnot(
  names(panels) == c("a_identity", "b_eigengene_matrix", "c_trajectories", "d_preservation")
)
plots <- purrr::map2(lapply(panels, `[[`, "panel"), letters[seq_along(panels)], tag_panel)
# Each row is wrapped so its panels line up with each other but not with the other row's axes.
figure <- wrap_plots(
  wrap_elements(full = wrap_plots(plots[1:2], widths = c(85, 95))),
  wrap_elements(full = wrap_plots(plots[3:4], widths = c(90, 90))),
  ncol = 1, heights = c(1, 1.05)
)
ggsave(
  here("04_Figures", "F5", "b_reports", "F5_composite.pdf"), figure,
  width = 180, height = 165, units = "mm", device = cairo_pdf
)

write_figure_workbook(panels, "F5")
