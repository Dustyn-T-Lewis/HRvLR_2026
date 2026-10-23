# S1 e: the T2 and T3 levels, the two PCAs F1 b leaves out, drawn the same way.
library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pca_grid.R"))

inputs <- pca_inputs
holds <- paste(
  "Panel e: PC1 and PC2 per person at T2 and T3, with each time point's PERMANOVA for group,",
  "gap k and forced-split adjusted Rand index."
)
panel_data <- pca_data(c("T2", "T3"))

panel <- pca_row(
  panel_data,
  title = "Proteome at T2 and T3",
  subtitle = paste(
    "As Figure 1b for the trained and post-bout levels; at both, the forced second cluster",
    "is S28's thin biopsy alone"
  )
)

save_panel(panel, "S1", "e_pca_levels", 125, 72)
