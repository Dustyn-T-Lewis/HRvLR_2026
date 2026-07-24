# HR/LR classification at the T1 contrast (pam engine), across feature spaces.
source(here::here("functions", "pred_leaf.R"))
run_class_leaf(
  pred_load(),
  contrast = "T1",
  model = "pam",
  out_dir = here::here(
    "04_Figures", "F05_prediction", "group_HRvLR", "T1", "pam"
  )
)
