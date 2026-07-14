#!/usr/bin/env Rscript
# Training concordance: HR vs LR over T1 to T2. Where the two groups adapt
# together and where they diverge. Rendered by the shared concordance driver.
pacman::p_load(here)
source(here("04_Figures", "F00_PILOT", "concordance_training", "a_script", "concordance_figure.R"))

render_concordance_figure(list(
  fig_id = "concordance_training",
  fig_dir = here("04_Figures", "F00_PILOT", "concordance_training"),
  c_hi = "Training_HR", c_lo = "Training_LR", c_int = "Training_Interaction",
  lo_levels = c("LR_T1", "LR_T2"),
  labels = list(
    hi = "HR", lo = "LR", phase = "Training",
    x = "HR training logFC", y = "LR training logFC",
    x_short = "HR training", y_short = "LR training"
  ),
  title = "Figure 4. Training response: shared adaptation and responder-specific divergence",
  subtitle = "HR vs LR over T1 to T2. Diagonal = concordant; off-diagonal and orange = interaction-significant divergence."
))
