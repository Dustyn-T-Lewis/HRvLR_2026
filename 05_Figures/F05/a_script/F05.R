# F05: classification and association at the protein, set and module levels, read against chance.
source(here::here("05_Figures", "theme.R"))
library(patchwork)

panels <- load_panels("F05")
composite <- panels$A_classification / panels$B_association +
  plot_layout(heights = c(1, 1.1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
save_composite(composite, panels, "F05", width = 178, height = 150)
