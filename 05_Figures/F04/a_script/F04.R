# F04: networks. Panels from the 04_Network workbooks.
source(here::here("05_Figures", "theme.R"))
library(patchwork)

panels <- load_panels("F04")
composite <- (panels$A_modules | panels$B_string) /
  (panels$C_preservation | panels$D_eigengene_tests) +
  plot_annotation(tag_levels = "A")
save_composite(composite, panels, "F04", width = 178, height = 140)
