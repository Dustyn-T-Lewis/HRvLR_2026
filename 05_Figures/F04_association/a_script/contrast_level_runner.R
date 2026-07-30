# One contrast heatmap per feature level. Each level owns its subdirectory, so a
# level can be rebuilt and read on its own while the composite still assembles
# from the same three panels.

pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "contrast_heatmap.R"))

CONTRAST_PANEL_HEIGHT <- c(modules = 132, pathways = 148, proteins = 148)

# The composite drops the caveat block, roughly 30 mm of text per level.
COMPOSITE_PANEL_HEIGHT <- c(modules = 104, pathways = 118, proteins = 118)
CONTRAST_PANEL_WIDTH <- 210

run_contrast_level <- function(level, caption = TRUE) {
  built <- build_contrast_level(level, caption = caption)
  stem <- here("05_Figures", "F04_association", level)
  for (sub in c("b_reports", "c_data")) {
    dir.create(file.path(stem, sub), recursive = TRUE, showWarnings = FALSE)
  }

  save_panel(
    built$panel, file.path(stem, "b_reports", paste0("F04_", level)),
    width = CONTRAST_PANEL_WIDTH, height = CONTRAST_PANEL_HEIGHT[[level]]
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
