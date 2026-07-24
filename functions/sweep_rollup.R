# Root roll-up: gather every leaf's metric summary into one raw table for later
# review. This is a data aggregation, not an interpretive composite -- it is the
# substrate the manual method-selection step reads. Leaves are the source of
# truth; the roll-up is rebuilt from them and never written concurrently with a
# leaf run.

pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "sweep_grid.R"))

rollup_root <- function(root, p_col) {
  files <- Sys.glob(
    file.path(sweep_root_dir(root), "*", "*", "*", "c_data", "results.xlsx")
  )
  all_cells <- bind_rows(lapply(files, function(f) {
    openxlsx::read.xlsx(f, sheet = "summary")
  }))
  leads <- all_cells |>
    filter(.data$B == max(.data$B), .data[[p_col]] < 0.05) |>
    arrange(.data[[p_col]])

  dir.create(file.path(sweep_root_dir(root), "c_data"),
    recursive = TRUE, showWarnings = FALSE
  )
  write_sweep_workbook(
    file.path(sweep_root_dir(root), "c_data", "results.xlsx"),
    list(all_cells = all_cells, leads = leads)
  )
  message(sprintf(
    "%s roll-up: %d leaves, %d cells, %d nominal leads (p<.05 at max B)",
    root, length(files), nrow(all_cells), nrow(leads)
  ))
  invisible(all_cells)
}
