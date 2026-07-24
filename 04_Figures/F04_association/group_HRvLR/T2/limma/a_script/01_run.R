# HR/LR differential abundance at the T2 contrast, limma engine, across spaces.
source(here::here("functions", "assoc_leaf.R"))
run_group_leaf(
  assoc_load(),
  contrast = "T2",
  method = "limma",
  out_dir = here::here("04_Figures", "F04_association", "group_HRvLR", "T2", "limma")
)
