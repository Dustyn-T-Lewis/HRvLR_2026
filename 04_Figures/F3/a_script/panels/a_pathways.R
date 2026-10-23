# Pathways per block: Training and Acute split by which group fgsea calls (shared, HR only, LR
# only), and the between-group window. Each tree cluster at Jaccard height 0.85
# (02_HRvLR/02_Set_Contrasts) shows its most specific set, the first `cap` per block by rank, drawn
# in every contrast with fgsea NES and FDR. The interaction window is in S3.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))

inputs <- c(set_tests = "02_HRvLR/02_Set_Contrasts/c_data/set_tests.rds")
cap <- 8
holds <- paste(
  sprintf("Per block, the most specific set of each tree cluster at height 0.85 (up to %d", cap),
  "by rank) in every contrast: fgsea NES and FDR within the database."
)
built <- cluster_grid(
  readRDS(here(inputs[["set_tests"]])), "0.85",
  cap = cap, windows = c("Training", "Acute", "Between groups"),
  title = "Which pathways each group moves, and where HR and LR differ",
  subtitle = sprintf(
    paste0(
      "Training and acute: fgsea FDR < %s in both groups, HR only or LR only (candidate sets ",
      "in brackets); between groups: two of three time points\none set per tree cluster ",
      "(Jaccard height 0.85), up to %d per block; * calls of opposite sign; ",
      "\u2020 the interaction also reaches FDR < %s"
    ),
    alpha_fdr, cap, alpha_fdr
  )
)
panel <- built$panel
panel_data <- built$panel_data
panel_height <- built$height

save_panel(panel, "F3", "a_pathways", 170, panel_height)
