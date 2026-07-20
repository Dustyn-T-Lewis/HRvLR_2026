# Leaf drivers for the prediction suite. Each feature-family leaf sources this
# and calls run_class_arm() or run_cont_arm() with its feature space and
# output directory, so the phase x model (x outcome) sweep, permutation
# nulls, BH correction, source workbook and panel all live in one place.

pacman::p_load(here, dplyr, readr, ggplot2, openxlsx)

source(here("functions", "shared_style.R"))
source(here("04_Figures", "F05_prediction", "shared", "a_script", "_features.R"))
source(here("functions", "shared_prediction.R"))
source(here("04_Figures", "F05_prediction", "shared", "a_script", "_panels.R"))

PRED_MODELS <- c("glmnet", "spls")
PRED_PHASES <- c("baseline", "training", "acute")
CONT_OUTCOMES <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)

leaf_dirs <- function(out_dir) {
  dat <- file.path(out_dir, "c_data")
  rpt <- file.path(out_dir, "b_reports", "panels")
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

# Responder-class (HR vs LR) prediction for one feature space, sweeping
# phase x model. One metrics table with a phase column, BH applied once.
run_class_arm <- function(bundle, out_dir, space, nperm = N_PERM,
                          phases = PRED_PHASES) {
  dirs <- leaf_dirs(out_dir)
  y_all <- pred_outcome(bundle, "group")

  cells <- expand.grid(
    phase = phases, model = PRED_MODELS,
    stringsAsFactors = FALSE
  )
  summ <- list()
  roc <- list()
  preds <- list()
  for (i in seq_len(nrow(cells))) {
    ph <- cells$phase[i]
    mo <- cells$model[i]
    al <- align_xy(
      pred_phase_matrix(bundle$feature_sets[[space]], bundle$meta, ph),
      y_all
    )
    message(sprintf(
      "[class %s] %s x %s  (n=%d, p=%d)",
      space, ph, mo, length(al$y), ncol(al$x)
    ))
    res <- run_class_cell(al$x, al$y, mo, nperm = nperm)
    summ[[i]] <- cbind(space = space, phase = ph, res$summary)
    roc[[i]] <- cbind(space = space, phase = ph, res$roc)
    preds[[i]] <- cbind(space = space, phase = ph, res$preds)
  }

  summ <- bind_rows(summ) |> mutate(q_bh = p.adjust(perm_p, "BH"))
  roc <- bind_rows(roc)
  preds <- bind_rows(preds)

  write_leaf_workbook(
    file.path(dirs$dat, "class_source_data.xlsx"),
    list(metrics = summ, roc = roc, predictions = preds)
  )

  panel <- build_class_panel(summ, roc, space)
  ggsave(file.path(dirs$rpt, "class_prediction.png"), panel,
    width = 190, height = 150, units = "mm", dpi = 300
  )
  ggsave(file.path(dirs$rpt, "class_prediction.pdf"), panel,
    width = 190, height = 150, units = "mm", device = PDF_DEVICE
  )

  message("\n=== Responder-class prediction (", space, ") ===")
  print(summ |> dplyr::select(phase, model, n, estimate, ci_lo, ci_hi, perm_p, q_bh))
  invisible(summ)
}

# Continuous-phenotype prediction for one feature space, sweeping
# outcome x phase x model. One metrics table with a phase column, BH
# applied once.
run_cont_arm <- function(bundle, out_dir, space, nperm = N_PERM,
                         phases = PRED_PHASES) {
  dirs <- leaf_dirs(out_dir)

  cells <- expand.grid(
    outcome = CONT_OUTCOMES, phase = phases,
    model = PRED_MODELS, stringsAsFactors = FALSE
  )
  summ <- list()
  preds <- list()
  for (i in seq_len(nrow(cells))) {
    oc <- cells$outcome[i]
    ph <- cells$phase[i]
    mo <- cells$model[i]
    y_all <- pred_outcome(bundle, oc)
    al <- align_xy(
      pred_phase_matrix(bundle$feature_sets[[space]], bundle$meta, ph),
      y_all
    )
    message(sprintf(
      "[cont %s] %s x %s x %s  (n=%d, p=%d)",
      space, oc, ph, mo, length(al$y), ncol(al$x)
    ))
    res <- run_cont_cell(al$x, al$y, mo, oc, nperm = nperm)
    summ[[i]] <- cbind(space = space, phase = ph, res$summary)
    preds[[i]] <- cbind(space = space, phase = ph, res$preds)
  }

  summ <- bind_rows(summ) |> mutate(q_bh = p.adjust(perm_p_q2, "BH"))
  preds <- bind_rows(preds)

  write_leaf_workbook(
    file.path(dirs$dat, "cont_source_data.xlsx"),
    list(metrics = summ, predictions = preds)
  )

  panel <- build_cont_panel(summ, preds, space)
  ggsave(file.path(dirs$rpt, "cont_prediction.png"), panel,
    width = 260, height = 230, units = "mm", dpi = 300
  )
  ggsave(file.path(dirs$rpt, "cont_prediction.pdf"), panel,
    width = 260, height = 230, units = "mm", device = PDF_DEVICE
  )

  message("\n=== Continuous-phenotype prediction (", space, ") ===")
  print(summ |> dplyr::select(
    outcome, phase, model, n, q2, rmse, spearman,
    perm_p_q2, q_bh
  ))
  invisible(summ)
}
