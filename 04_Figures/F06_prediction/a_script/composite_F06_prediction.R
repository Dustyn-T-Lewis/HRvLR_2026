suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
}))

panel <- build_cont_composite("F06_prediction")
dir.create(here("04_Figures", "F06_prediction", "b_reports"),
  recursive = TRUE, showWarnings = FALSE
)
save_panel(
  panel,
  here("04_Figures", "F06_prediction", "b_reports", "F06_prediction_composite"),
  width = 250, height = 320
)
