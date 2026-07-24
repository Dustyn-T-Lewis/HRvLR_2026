# Leaf driver for the F05 prediction suite. Each terminal leaf is one contrast x
# one model; it runs that cell across the three feature spaces (panels) and
# the source workbook, and renders the leaf panel. BH is deferred to composite
# assembly because it is applied once across the whole figure grid, so each leaf
# stores the raw permutation p and the composite adds the BH q.

pacman::p_load(here, dplyr, openxlsx, ggplot2)
source(here("functions", "shared_style.R"))
source(here("functions", "pred_features.R"))
source(here("functions", "shared_prediction.R"))
source(here("functions", "pred_panels.R"))

PRED_SPACES <- c("proteins", "singscore", "eigengenes")

# Reading label for each classification contrast (methods doc Part 0, Part 8).
CONTRAST_ROLE <- c(
  T1 = "prospective forecast", T2 = "separability", T3 = "separability",
  training = "trajectory (T2-T1)", acute = "trajectory (T3-T2)",
  baseline = "prospective forecast"
)

CONT_OUTCOMES <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)

leaf_dirs <- function(out_dir) {
  dat <- file.path(out_dir, "c_data")
  rpt <- file.path(out_dir, "b_reports")
  dir.create(dat, recursive = TRUE, showWarnings = FALSE)
  dir.create(rpt, recursive = TRUE, showWarnings = FALSE)
  list(dat = dat, rpt = rpt)
}

write_leaf_workbook <- function(path, sheets) {
  wb <- createWorkbook()
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]])
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}

tag_cell <- function(df, contrast, space) {
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }
  cbind(contrast = contrast, space = space, df)
}

# Classification leaf: contrast x model across the three feature spaces.
run_class_leaf <- function(bundle, contrast, model, out_dir,
                           nperm = N_PERM, spaces = PRED_SPACES) {
  dirs <- leaf_dirs(out_dir)
  y_all <- pred_outcome(bundle, "group")

  cells <- lapply(spaces, function(sp) {
    al <- align_xy(
      pred_contrast_matrix(bundle$feature_sets[[sp]], bundle$meta, contrast),
      y_all
    )
    message(sprintf(
      "[class] %s x %s x %s  (n=%d, p=%d)",
      contrast, model, sp, length(al$y), ncol(al$x)
    ))
    res <- run_class_cell(al$x, al$y, model, nperm = nperm)
    list(
      summary = tag_cell(res$summary, contrast, sp),
      roc = tag_cell(res$roc, contrast, sp),
      preds = tag_cell(res$preds, contrast, sp),
      selection = tag_cell(res$selection, contrast, sp)
    )
  })

  summ <- bind_rows(lapply(cells, `[[`, "summary"))
  roc <- bind_rows(lapply(cells, `[[`, "roc"))
  preds <- bind_rows(lapply(cells, `[[`, "preds"))
  selection <- bind_rows(lapply(cells, `[[`, "selection"))

  write_leaf_workbook(
    file.path(dirs$dat, "source_data.xlsx"),
    list(metrics = summ, roc = roc, predictions = preds, selection = selection)
  )

  panel <- build_class_leaf_panel(summ, roc, selection, contrast, model)
  save_panel(panel, file.path(dirs$rpt, "panel"), width = 190, height = 130)
  invisible(summ)
}

# Continuous leaf: model across the three feature spaces, sweeping the six
# adaptation outcomes at the baseline phase.
run_cont_leaf <- function(bundle, model, out_dir, nperm = N_PERM,
                          spaces = PRED_SPACES, outcomes = CONT_OUTCOMES,
                          contrast = "baseline") {
  dirs <- leaf_dirs(out_dir)

  cells <- list()
  for (sp in spaces) {
    fmat <- pred_contrast_matrix(
      bundle$feature_sets[[sp]], bundle$meta, contrast
    )
    for (oc in outcomes) {
      al <- align_xy(fmat, pred_outcome(bundle, oc))
      message(sprintf(
        "[cont] %s x %s x %s  (n=%d, p=%d)",
        oc, model, sp, length(al$y), ncol(al$x)
      ))
      res <- run_cont_cell(al$x, al$y, model, oc, nperm = nperm)
      cells[[length(cells) + 1]] <- list(
        summary = cbind(space = sp, res$summary),
        preds = cbind(space = sp, res$preds)
      )
    }
  }

  summ <- bind_rows(lapply(cells, `[[`, "summary"))
  preds <- bind_rows(lapply(cells, `[[`, "preds"))

  write_leaf_workbook(
    file.path(dirs$dat, "source_data.xlsx"),
    list(metrics = summ, predictions = preds)
  )

  panel <- build_cont_leaf_panel(summ, preds, model)
  save_panel(panel, file.path(dirs$rpt, "panel"), width = 230, height = 150)
  invisible(summ)
}
