# Pathways per window: the representative of each tree cluster at Jaccard height 0.85
# (02_HRvLR/02_Set_Contrasts), up to ten per window by rank, drawn in all nine contrasts with fgsea
# NES and FDR.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))

inputs <- c(set_tests = "02_HRvLR/02_Set_Contrasts/c_data/set_tests.rds")
holds <- paste(
  "Per window, the representative set of each tree cluster at height 0.85 (up to ten by rank)",
  "in all nine contrasts: fgsea NES and FDR within the database."
)
built <- cluster_grid(
  readRDS(here(inputs[["set_tests"]])), "0.85",
  cap = 10,
  title = "Pathways per window, one per cluster of overlapping sets",
  subtitle = paste0(
    "fgsea FDR < 0.05 in both groups of a window (between groups: two of three time points; ",
    "interaction: either)\none set per tree cluster (Jaccard height 0.85), up to ten per window, ",
    "shown in all nine contrasts; * calls of opposite sign"
  )
)
panel <- built$panel
panel_data <- built$panel_data
panel_height <- built$height

ggsave(
  here("04_Figures", "F3", "b_reports", "panels", "a_pathways.pdf"), panel,
  width = 180, height = panel_height, units = "mm", device = cairo_pdf
)
