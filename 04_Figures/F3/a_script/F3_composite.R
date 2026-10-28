# F3: the pathways behind F1 e's counts, one representative per cluster of overlapping sets.
library(here)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F3")
stopifnot(names(panels) == "a_pathways")
height <- panels$a_pathways$panel_height
stopifnot(height <= 240)
save_figure(tag_panel(panels$a_pathways$panel, "a"), "F3", height)

write_figure_workbook(panels, "F3")
