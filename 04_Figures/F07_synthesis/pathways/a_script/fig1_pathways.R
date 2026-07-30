# Rebuild the pathways panel of Figure 1 on its own, without the other levels.

suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "synthesis_heatmap.R"))
  source(here(
    "04_Figures", "F07_synthesis", "a_script", "fig1_level_runner.R"
  ))
}))

invisible(run_fig1_level("pathways"))
