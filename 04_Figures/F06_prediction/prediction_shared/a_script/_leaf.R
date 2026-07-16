# Leaf drivers for the prediction suite. Each phase leaf sources this and calls
# run_class_leaf() or run_cont_leaf() with its phase and output directory, so
# the feature-space x model x outcome sweep, permutation nulls, BH correction,
# metric CSVs and panels all live in one place.

pacman::p_load(here, dplyr, readr, ggplot2)

source(here("04_Figures", "functions", "shared_style.R"))
source(here("04_Figures", "F06_prediction", "prediction_shared", "a_script", "_features.R"))
source(here("04_Figures", "F06_prediction", "prediction_shared", "a_script", "_harness.R"))
source(here("04_Figures", "F06_prediction", "prediction_shared", "a_script", "_panels.R"))

PRED_SPACES <- c("singscore", "proteins", "eigengenes")
PRED_MODELS <- c("glmnet", "spls")
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

run_class_leaf <- function(bundle, phase, out_dir, nperm = N_PERM) {
  dirs <- leaf_dirs(out_dir)
  y_all <- pred_outcome(bundle, "group")

  cells <- expand.grid(
    space = PRED_SPACES, model = PRED_MODELS,
    stringsAsFactors = FALSE
  )
  summ <- list()
  roc <- list()
  preds <- list()
  for (i in seq_len(nrow(cells))) {
    sp <- cells$space[i]
    mo <- cells$model[i]
    al <- align_xy(
      pred_phase_matrix(bundle$feature_sets[[sp]], bundle$meta, phase),
      y_all
    )
    message(sprintf(
      "[class %s] %s x %s  (n=%d, p=%d)",
      phase, sp, mo, length(al$y), ncol(al$x)
    ))
    res <- run_class_cell(al$x, al$y, mo, nperm = nperm)
    summ[[i]] <- cbind(phase = phase, space = sp, res$summary)
    roc[[i]] <- cbind(phase = phase, space = sp, res$roc)
    preds[[i]] <- cbind(phase = phase, space = sp, res$preds)
  }

  summ <- bind_rows(summ) |> mutate(q_bh = p.adjust(perm_p, "BH"))
  roc <- bind_rows(roc)
  preds <- bind_rows(preds)

  write_csv(summ, file.path(dirs$dat, "class_metrics.csv"))
  write_csv(roc, file.path(dirs$dat, "class_roc.csv"))
  write_csv(preds, file.path(dirs$dat, "class_predictions.csv"))

  panel <- build_class_panel(summ, roc, phase)
  ggsave(file.path(dirs$rpt, "class_prediction.png"), panel,
    width = 190, height = 150, units = "mm", dpi = 300
  )
  ggsave(file.path(dirs$rpt, "class_prediction.pdf"), panel,
    width = 190, height = 150, units = "mm", device = PDF_DEVICE
  )

  message("\n=== Responder-class prediction (", phase, ") ===")
  print(summ |> dplyr::select(space, model, n, estimate, ci_lo, ci_hi, perm_p, q_bh))
  invisible(summ)
}

run_cont_leaf <- function(bundle, phase, out_dir, nperm = N_PERM) {
  dirs <- leaf_dirs(out_dir)

  cells <- expand.grid(
    outcome = CONT_OUTCOMES, space = PRED_SPACES,
    model = PRED_MODELS, stringsAsFactors = FALSE
  )
  summ <- list()
  preds <- list()
  for (i in seq_len(nrow(cells))) {
    oc <- cells$outcome[i]
    sp <- cells$space[i]
    mo <- cells$model[i]
    y_all <- pred_outcome(bundle, oc)
    al <- align_xy(
      pred_phase_matrix(bundle$feature_sets[[sp]], bundle$meta, phase),
      y_all
    )
    message(sprintf(
      "[cont %s] %s x %s x %s  (n=%d, p=%d)",
      phase, oc, sp, mo, length(al$y), ncol(al$x)
    ))
    res <- run_cont_cell(al$x, al$y, mo, oc, nperm = nperm)
    summ[[i]] <- cbind(phase = phase, space = sp, res$summary)
    preds[[i]] <- cbind(phase = phase, space = sp, res$preds)
  }

  summ <- bind_rows(summ) |> mutate(q_bh = p.adjust(perm_p_q2, "BH"))
  preds <- bind_rows(preds)

  write_csv(summ, file.path(dirs$dat, "cont_metrics.csv"))
  write_csv(preds, file.path(dirs$dat, "cont_predictions.csv"))

  panel <- build_cont_panel(summ, preds, phase)
  ggsave(file.path(dirs$rpt, "cont_prediction.png"), panel,
    width = 260, height = 150, units = "mm", dpi = 300
  )
  ggsave(file.path(dirs$rpt, "cont_prediction.pdf"), panel,
    width = 260, height = 150, units = "mm", device = PDF_DEVICE
  )

  message("\n=== Continuous-phenotype prediction (", phase, ") ===")
  print(summ |> dplyr::select(
    outcome, space, model, n, q2, rmse, spearman,
    perm_p_q2, q_bh
  ))
  invisible(summ)
}
