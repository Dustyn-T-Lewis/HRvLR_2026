# HR/LR classification at the training contrast (glmnet engine), across feature spaces.
source(here::here("functions", "pred_leaf.R"))
run_class_leaf(
  pred_load(),
  contrast = "training",
  model = "glmnet",
  out_dir = here::here(
    "04_Figures", "F05_prediction", "group_HRvLR", "training", "glmnet"
  )
)
