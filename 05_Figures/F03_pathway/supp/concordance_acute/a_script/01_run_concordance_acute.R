#!/usr/bin/env Rscript
# Acute concordance: HR vs LR over T2 to T3. Same builders as the training arm on the
# acute contrast pair. Rendered by the shared concordance driver.
pacman::p_load(here)
source(here("05_Figures", "functions", "f00_concordance.R"))

render_concordance_figure(list(
  fig_id = "concordance_acute",
  fig_dir = here("05_Figures", "F03_pathway", "supp", "concordance_acute"),
  c_hi = "Acute_HR", c_lo = "Acute_LR", c_int = "Acute_Interaction",
  lo_levels = c("LR_T2", "LR_T3"),
  labels = list(
    hi = "HR", lo = "LR", phase = "Acute",
    x = "HR acute logFC", y = "LR acute logFC",
    x_short = "HR acute", y_short = "LR acute"
  ),
  title = "Pilot. Acute response: shared signal and responder-specific divergence",
  subtitle = "HR vs LR over T2 to T3. Diagonal = concordant; off-diagonal and orange = interaction-significant divergence."
))
