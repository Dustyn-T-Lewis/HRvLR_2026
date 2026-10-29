# S3 b: one set per tree cluster at Jaccard height 0.75 (02_HRvLR/02_Set_Contrasts), every
# cluster of every window; Figure 3 cuts at 0.85.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))

inputs <- c(set_tests = "02_HRvLR/02_Set_Contrasts/c_data/set_contrasts.rds")
list2env(cut_panel(readRDS(here(inputs[["set_tests"]])), "0.75"), environment())
save_panel(panel, "S3", "b_cut_075", 250, panel_height)
