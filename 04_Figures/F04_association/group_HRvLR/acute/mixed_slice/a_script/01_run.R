# HR/LR differential abundance at the acute contrast, mixed_slice engine, across spaces.
source(here::here("functions", "assoc_leaf.R"))
run_group_leaf(
  assoc_load(),
  contrast = "acute",
  method = "mixed_slice",
  out_dir = here::here("04_Figures", "F04_association", "group_HRvLR", "acute", "mixed_slice")
)
