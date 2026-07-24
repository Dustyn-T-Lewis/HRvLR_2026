# HR/LR classification at the T1 contrast (spls engine), across feature spaces.
source(here::here("functions", "pred_leaf.R"))
run_class_leaf(
  pred_load(),
  contrast = "T1",
  model = "spls",
  out_dir = here::here(
    "04_Figures", "F05_prediction", "group_HRvLR", "T1", "splsda"
  )
)
