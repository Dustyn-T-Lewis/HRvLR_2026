# F01: cohort and filtering. Panels from 00_Input/phenotype.csv and 01_Filtering's workbook.
source(here::here("05_Figures", "theme.R"))
library(patchwork)

panels <- load_panels("F01")
composite <- panels$A_phenotype / (panels$B_filtering | panels$C_blood) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "A")
save_composite(composite, panels, "F01", width = 178, height = 130)
