# One set per tree cluster at Jaccard height 0.75 (02_HRvLR/02_Set_Contrasts), every cluster of
# every block, drawn in all nine contrasts; Figure 3 draws the cut at 0.85.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))

inputs <- c(set_tests = "02_HRvLR/02_Set_Contrasts/c_data/set_tests.rds")
holds <- paste(
  "One set per tree cluster at height 0.75, every cluster, in all nine contrasts:",
  "fgsea NES and FDR within the database."
)
built <- cluster_grid(
  readRDS(here(inputs[["set_tests"]])), "0.75",
  title = "One set per cluster at tree height 0.75",
  subtitle = paste(
    "As Figure 3 with the tree cut at 0.75 instead of 0.85, every cluster shown;",
    "* calls of opposite sign; \u2020 interaction FDR < 0.05"
  )
)
panel <- built$panel
panel_data <- built$panel_data
panel_height <- built$height

ggsave(
  here("04_Figures", "S3", "b_reports", "panels", "b_cut_075.pdf"), panel,
  width = 180, height = panel_height, units = "mm", device = cairo_pdf
)
