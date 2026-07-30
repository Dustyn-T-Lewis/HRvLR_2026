#!/usr/bin/env Rscript
pacman::p_load(here, openxlsx)
source(here("functions", "sweep_split.R"))
source(here("functions", "sweep_cell_panel.R"))

root <- sweep_root_dir("F04_association")
leaves <- Sys.glob(file.path(root, "*", "*", "*", "c_data", "results.xlsx"))

for (path in leaves) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  n <- length(parts)
  method <- parts[[n - 2]]
  config <- parts[[n - 3]]
  level <- parts[[n - 4]]

  for (oc in split_assoc_leaf(path, root, level, config, method)) {
    d <- split_leaf_dir(root, level, config, oc, method)
    cell <- openxlsx::read.xlsx(file.path(d, "c_data", "results.xlsx"), "cell")
    panel <- build_assoc_cell_panel(cell, level, config, oc, method)
    save_panel(panel, file.path(d, "b_reports", "panel"),
      width = 210, height = 90
    )
  }
  message(sprintf("split %s/%s/%s", level, config, method))
}
