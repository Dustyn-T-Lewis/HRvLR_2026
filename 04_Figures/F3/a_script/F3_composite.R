library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F3")
stopifnot(names(panels) == "a_pathways")

height <- panels$a_pathways$panel_height
stopifnot(height <= 240)
figure <- panels$a_pathways$panel + labs(tag = "a") + figure_tags
ggsave(
  here("04_Figures", "F3", "b_reports", "F3_composite.pdf"), figure,
  width = 180, height = height, units = "mm", device = cairo_pdf
)

write_figure_workbook(panels, "F3")
