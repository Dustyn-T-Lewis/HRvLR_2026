# F5: modules. Row 1 reads each module across: name and size, STRING enrichment, eigengene tests.
# Row 2: the clearest module person by person, and preservation between the groups.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F5")
stopifnot(
  names(panels) == c("a_identity", "b_eigengene_matrix", "c_trajectories", "d_preservation")
)
identity <- panels$a_identity
matrix <- panels$b_eigengene_matrix$panel + ggplot2::theme(axis.text.y = ggplot2::element_blank())
stopifnot(identical(
  levels(identity$d$module), rev(panels$b_eigengene_matrix$labels$module)
))

# Widths on a 250 mm page. Row 1: size bars 100 mm (names included), STRING 35, eigengenes 115,
# one set of rows. Row 2: c 100 mm, d 150.
row_a_b <- wrap_plots(
  tag_panel(identity$size_plot, "a"), identity$string_plot, tag_panel(matrix, "b"),
  widths = c(100, 35, 115)
)
row_c_d <- wrap_plots(
  tag_panel(panels$c_trajectories$panel, "c"), tag_panel(panels$d_preservation$panel, "d"),
  widths = c(100, 150)
)
figure <- wrap_plots(
  wrap_elements(full = row_a_b), wrap_elements(full = row_c_d),
  ncol = 1, heights = c(80, 85)
)
save_figure(figure, "F5", 165)

write_figure_workbook(panels, "F5")
