# F2: the proteins each group changes with training or the bout. The exploratory HR - LR candidates
# are in S2.
library(here)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F2")
stopifnot(names(panels) == "a_responders")
figure <- tag_panel(panels$a_responders$panel, "a")
save_figure(figure, "F2", 140)

write_figure_workbook(panels, "F2")
