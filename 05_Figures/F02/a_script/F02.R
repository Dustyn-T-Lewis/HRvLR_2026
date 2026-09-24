# F02: the proteome. Panels from 02_Normalization's and 02_Differential's workbooks.
source(here::here("05_Figures", "theme.R"))
library(patchwork)

panels <- load_panels("F02")
composite <- (panels$A_pca | panels$B_hits) / (panels$C_p_histogram | panels$D_volcano) +
  plot_annotation(tag_levels = "A")
save_composite(composite, panels, "F02", width = 178, height = 140)
