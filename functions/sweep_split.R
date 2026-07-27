# Re-slice a swept leaf into one leaf per phenotype. The swept workbooks are the
# source of truth; nothing here recomputes a model.

pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "sweep_grid.R"))

PRED_SHEETS <- c("summary", "null", "predictions", "selection")

split_leaf_dir <- function(root, level, config, phenotype, model) {
  d <- file.path(root, level, config, phenotype, model)
  for (sub in c("b_reports", "c_data")) {
    dir.create(file.path(d, sub), recursive = TRUE, showWarnings = FALSE)
  }
  d
}

read_sheet <- function(path, sheet) {
  openxlsx::read.xlsx(path, sheet)
}

# A sheet is filtered only when it carries the outcome key; F05's null,
# predictions and selection sheets do not, so they copy across whole.
filter_to_outcome <- function(df, outcome) {
  if (is.null(df) || !"outcome" %in% names(df)) {
    return(df)
  }
  df[df$outcome == outcome, , drop = FALSE]
}

split_pred_leaf <- function(path, root, level, config, method,
                            phenotype = NULL) {
  sheets <- lapply(PRED_SHEETS, function(s) {
    if (s %in% openxlsx::getSheetNames(path)) read_sheet(path, s)
  })
  names(sheets) <- PRED_SHEETS

  outcomes <- if (!is.null(phenotype)) {
    phenotype
  } else if ("outcome" %in% names(sheets$summary)) {
    unique(sheets$summary$outcome)
  } else {
    stop("no outcome column and no phenotype given for ", path)
  }

  for (oc in outcomes) {
    keep <- if (is.null(phenotype)) {
      lapply(sheets, filter_to_outcome, outcome = oc)
    } else {
      sheets
    }
    d <- split_leaf_dir(root, level, config, oc, method)
    write_sweep_workbook(
      file.path(d, "c_data", "results.xlsx"),
      Filter(Negate(is.null), keep)
    )
  }
  outcomes
}

split_assoc_leaf <- function(path, root, level, config, method) {
  present <- openxlsx::getSheetNames(path)
  outcomes <- setdiff(present, "summary")
  summary_all <- read_sheet(path, "summary")

  for (oc in outcomes) {
    d <- split_leaf_dir(root, level, config, oc, method)
    write_sweep_workbook(
      file.path(d, "c_data", "results.xlsx"),
      list(
        cell = read_sheet(path, oc),
        summary = summary_all[summary_all$outcome == oc, , drop = FALSE]
      )
    )
  }
  outcomes
}
