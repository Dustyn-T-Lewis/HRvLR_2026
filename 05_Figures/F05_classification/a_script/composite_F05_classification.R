suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
}))

panel <- build_class_composite("F05_classification")
dir.create(here("05_Figures", "F05_classification", "b_reports"),
  recursive = TRUE, showWarnings = FALSE
)
save_panel(
  panel,
  here(
    "05_Figures", "F05_classification", "b_reports",
    "F05_classification_composite"
  ),
  width = 250, height = 300
)
