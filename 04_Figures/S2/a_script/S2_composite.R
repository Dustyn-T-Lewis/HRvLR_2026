library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S2")
stopifnot(names(panels) == "a_auc_top")
figure <- panels$a_auc_top$panel + labs(tag = "a") + figure_tags
ggsave(
  here("04_Figures", "S2", "b_reports", "S2_composite.pdf"), figure,
  width = 180, height = 140, units = "mm", device = cairo_pdf
)

write_figure_workbook(panels, "S2")
