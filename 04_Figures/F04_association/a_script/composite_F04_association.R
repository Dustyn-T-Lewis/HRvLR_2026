suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
}))

panel <- build_assoc_composite("F04_association")
dir.create(here("04_Figures", "F04_association", "b_reports"),
  recursive = TRUE, showWarnings = FALSE
)
save_panel(
  panel,
  here("04_Figures", "F04_association", "b_reports", "F04_association_composite"),
  width = 320, height = 210
)
