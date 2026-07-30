# All three feature levels on one sheet. patchwork drops a nested plot's own
# subtitle, so each level is wrapped before stacking.
#
# The composite carries titles and subtitles only. The caveat block is identical
# across the three levels, and stacking three copies of it buries the panels it
# belongs to. It stays on the per-level figures, the ones read closely.

pacman::p_load(here, patchwork, ggplot2, openxlsx, dplyr)
source(here(
  "05_Figures", "F04_association", "a_script", "contrast_level_runner.R"
))

LEVEL_ORDER <- c("proteins", "modules", "pathways")

build_f04_composite <- function() {
  built <- lapply(LEVEL_ORDER, run_contrast_level)
  names(built) <- LEVEL_ORDER
  bare <- lapply(LEVEL_ORDER, build_contrast_level, caption = FALSE)
  heights <- COMPOSITE_PANEL_HEIGHT[LEVEL_ORDER]
  panel <- Reduce(`/`, lapply(bare, function(b) wrap_elements(b$panel))) +
    plot_layout(heights = unname(heights))
  list(panel = panel, built = built, height = sum(heights))
}

composite <- build_f04_composite()
stem <- here("05_Figures", "F04_association")

save_panel(
  composite$panel, file.path(stem, "b_reports", "MAIN_F04_association"),
  width = CONTRAST_PANEL_WIDTH, height = composite$height
)

# One workbook for the whole stage: every contrast at every level, plus the row
# sets each panel actually shows.
write.xlsx(
  c(
    lapply(composite$built, function(b) b$data),
    stats::setNames(
      lapply(composite$built, function(b) select(b$rows, -"fill")),
      paste0(LEVEL_ORDER, "_rows")
    )
  ),
  file.path(stem, "c_data", "F04_association_source_data.xlsx")
)
