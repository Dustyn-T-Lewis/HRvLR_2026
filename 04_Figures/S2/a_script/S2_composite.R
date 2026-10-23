# S2: single-feature classification, one panel across the page.
library(here)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S2")
stopifnot(names(panels) == "a_auc_top")
figure <- tag_panel(panels$a_auc_top$panel, "a")
save_figure(figure, "S2", 120)

write_figure_workbook(panels, "S2")
