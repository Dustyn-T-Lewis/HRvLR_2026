#!/usr/bin/env Rscript
# F05 Acute concordance: HR vs LR over T2 to T3. Same builders as F04 on the
# acute contrast pair. Rendered by the shared concordance driver.
pacman::p_load(here)
source(here("04_Figures", "functions", "concordance_figure.R"))

render_concordance_figure(list(
  fig_id = "F05",
  fig_dir = here("04_Figures", "extras", "concordance_acute"),
  c_hi = "Acute_HR", c_lo = "Acute_LR", c_int = "Acute_Interaction",
  lo_levels = c("LR_T2", "LR_T3"),
  labels = list(
    hi = "HR", lo = "LR", phase = "Acute",
    x = "HR acute logFC", y = "LR acute logFC",
    x_short = "HR acute", y_short = "LR acute"
  ),
  title = "Figure 5. Acute response: shared signal and responder-specific divergence",
  subtitle = "HR vs LR over T2 to T3. Diagonal = concordant; off-diagonal and orange = interaction-significant divergence."
))
