# Figure 1: feature x trait heatmaps, one per feature level, plus the stacked
# composite. Each level writes its panel and its tile table into its own
# subdirectory so the rows behind a figure stay auditable next to it.

suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "synthesis_heatmap.R"))
  source(here(
    "04_Figures", "F07_synthesis", "a_script", "fig1_level_runner.R"
  ))
  pacman::p_load(dplyr, patchwork)
}))

built <- lapply(SYNTH_LEVELS, run_fig1_level)
names(built) <- SYNTH_LEVELS

# Each level is wrapped before stacking. Nesting a patchwork that carries its
# own annotations drops them, which would strip the survivor counts, the
# correction statement and the row-name footnote from the composite.
composite <- wrap_plots(
  lapply(built, function(b) wrap_elements(b$panel)),
  ncol = 1
) +
  plot_layout(heights = unname(FIG1_HEIGHT[SYNTH_LEVELS])) +
  plot_annotation(
    title = "Figure 1  Feature x adaptation associations at three levels",
    subtitle = paste(
      "Spearman rho against six adaptation outcomes, nested in seven timepoint",
      "configs. Value printed in tile, stars mark nominal p, black outline",
      "marks BH q < .05 applied within the cell.\nOne method throughout:",
      "Spearman is fixed on design grounds, not chosen on yield. The",
      "per-method consequence is a separate panel."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 4),
      plot.subtitle = element_text(
        face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
      )
    )
  )

save_panel(
  composite,
  file.path(sweep_root_dir(FIG1_STAGE), "b_reports", "FIG1_heatmap_composite"),
  340, sum(FIG1_HEIGHT) + 20
)
message("composite written to ", file.path(
  sweep_root_dir(FIG1_STAGE), "b_reports"
))
