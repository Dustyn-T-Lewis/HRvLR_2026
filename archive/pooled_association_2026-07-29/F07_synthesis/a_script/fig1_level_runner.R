# Writes one level's Figure 1 panel and its tile table into that level's own
# subdirectory. The stage driver loops this; each level's a_script calls it for
# its own level so a single level can be rebuilt without the other two.

pacman::p_load(here, dplyr, openxlsx)

FIG1_STAGE <- "F07_synthesis"
FIG1_HEIGHT <- c(modules = 118, pathways = 150, proteins = 150)

fig1_level_dir <- function(level, sub) {
  d <- file.path(sweep_root_dir(FIG1_STAGE), level, sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

run_fig1_level <- function(level) {
  out <- build_fig1_level(level)
  save_panel(
    out$panel, file.path(fig1_level_dir(level, "b_reports"), "fig1_heatmap"),
    340, FIG1_HEIGHT[[level]]
  )
  write_sweep_workbook(
    file.path(fig1_level_dir(level, "c_data"), "fig1_heatmap.xlsx"),
    list(
      tiles = arrange(out$data, .data$p),
      rows = mutate(out$rows, feature = as.character(.data$feature))
    )
  )
  message(level, ": ", nrow(out$rows), " rows, ", nrow(out$data), " tiles")
  out
}
