# Root roll-up: which cells lead, read off the store. The summary table already
# holds every cell, so this no longer gathers or writes anything -- it derives
# the lead set and reports the screen it came from.

pacman::p_load(here, dplyr)
source(here("functions", "sweep_grid.R"))

root_leads <- function(root, p_col) {
  kind <- root_kind(root)
  root_cells(root) |>
    filter(is_lead(.data[[root_metric_col(kind)]], .data[[p_col]], kind)) |>
    arrange(.data[[p_col]])
}

rollup_root <- function(root, p_col) {
  all_cells <- read_sweep_store(root, "summary")
  leads <- root_leads(root, p_col)
  message(sprintf(
    "%s roll-up: %d cells, %d nominal leads (p<.05 at max B)",
    root, nrow(all_cells), nrow(leads)
  ))
  invisible(all_cells)
}
