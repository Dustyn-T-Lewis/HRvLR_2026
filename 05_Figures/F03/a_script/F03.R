# F03: pathways. Panels from 01_run_fgsea_and_fry's and 03_enrich_scatter_fgsea's workbooks.
source(here::here("05_Figures", "theme.R"))
library(patchwork)

panels <- load_panels("F03")
composite <- (panels$A_set_counts | panels$B_nes_training) / panels$C_top_sets +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A")
save_composite(composite, panels, "F03", width = 178, height = 135)
