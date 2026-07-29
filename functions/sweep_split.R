# Re-slice a swept leaf into one leaf per phenotype. The swept workbooks are the
# source of truth; nothing here recomputes a model.

pacman::p_load(here, openxlsx)
source(here("functions", "sweep_grid.R"))

PRED_SHEETS <- c("summary", "null", "predictions", "selection")

# `root` is an already-resolved figure directory, so callers pass
# sweep_root_dir("F06_prediction") rather than the bare name.
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

# F04's per-outcome counts sheet is `cell_summary`, not `summary`: `summary`
# means the metric row in F05 and F06, and one name cannot mean both without
# silently handing a generic reader the wrong table.
split_assoc_leaf <- function(path, root, level, config, method) {
  present <- openxlsx::getSheetNames(path)
  outcomes <- setdiff(present, "cell_summary")
  summary_all <- read_sheet(path, "cell_summary")

  for (oc in outcomes) {
    d <- split_leaf_dir(root, level, config, oc, method)
    write_sweep_workbook(
      file.path(d, "c_data", "results.xlsx"),
      list(
        cell = read_sheet(path, oc),
        cell_summary = summary_all[summary_all$outcome == oc, , drop = FALSE]
      )
    )
  }
  outcomes
}
