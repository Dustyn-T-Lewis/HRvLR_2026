#!/usr/bin/env Rscript
pacman::p_load(here, openxlsx)
source(here("functions", "sweep_split.R"))
source(here("functions", "sweep_cell_panel.R"))

root <- sweep_root_dir("F06_prediction")
leaves <- Sys.glob(file.path(root, "*", "*", "*", "c_data", "results.xlsx"))

for (path in leaves) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  n <- length(parts)
  method <- parts[[n - 2]]
  config <- parts[[n - 3]]
  level <- parts[[n - 4]]

  for (oc in split_pred_leaf(path, root, level, config, method)) {
    d <- split_leaf_dir(root, level, config, oc, method)
    out <- file.path(d, "c_data", "results.xlsx")
    sheets <- lapply(PRED_SHEETS, function(s) openxlsx::read.xlsx(out, s))
    names(sheets) <- PRED_SHEETS
    panel <- build_cont_cell_panel(sheets, level, config, oc, method)
    save_panel(panel, file.path(d, "b_reports", "panel"),
      width = 220, height = 110
    )
  }
  message(sprintf("split %s/%s/%s", level, config, method))
}
