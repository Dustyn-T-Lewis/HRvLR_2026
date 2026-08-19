# One contrast heatmap per feature level. Each level owns its subdirectory, so a
# level can be rebuilt and read on its own while the composite still assembles
# from the same three panels.

pacman::p_load(here, dplyr, openxlsx)
source(here("04_Figures", "F04_association", "a_script", "contrast_heatmap.R"))

COMPOSITE_PANEL_HEIGHT <- c(modules = 104, pathways = 118, proteins = 118)
CONTRAST_PANEL_WIDTH <- 210

# The caveats ship as a legend file rather than baked into the render. A journal
# wants the legend as text, and 30 mm of 5.4 pt grey under every panel buries
# the thing it explains.
write_legend <- function(built, level, path) {
  writeLines(
    c(
      sprintf(
        "**%s: how high responders differ from low responders.**",
        str_to_sentence(level)
      ),
      "",
      contrast_subtitle(built$data, built$rows, level),
      "",
      contrast_caption(level, built$rows)
    ),
    path
  )
}

run_contrast_level <- function(level) {
  built <- build_contrast_level(level, caption = FALSE)
  stem <- here("04_Figures", "F04_association", level)
  for (sub in c("b_reports", "c_data")) {
    dir.create(file.path(stem, sub), recursive = TRUE, showWarnings = FALSE)
  }

  save_panel(
    built$panel, file.path(stem, "b_reports", paste0("F04_", level)),
    width = CONTRAST_PANEL_WIDTH, height = COMPOSITE_PANEL_HEIGHT[[level]]
  )
  write_legend(
    built, level,
    file.path(stem, "b_reports", paste0("F04_", level, "_legend.md"))
  )
  write.xlsx(
    list(
      shown = filter(
        built$data, .data$feature %in% levels(built$rows$feature)
      ),
      rows = select(built$rows, -"fill"),
      all_contrasts = built$data
    ),
    file.path(stem, "c_data", paste0("F04_", level, "_source_data.xlsx"))
  )
  built
}
