# Leakage-free nested leave-one-subject-out harness for F04_classification
# and F04_association. Rows are subjects, so leave-one-subject-out is
# leave-one-row-out. Every fold: features are z-scored on the training
# subjects and that centre/scale is applied to the held-out subject; the
# elastic-net alpha and lambda are tuned by an inner LOSO on the training
# subjects only; the held-out subject is never seen until it is predicted.
# Every fit also reports which features it selected, so the outer folds
# yield a per-feature selection frequency. The permutation null shuffles
# the outcome across subjects and re-runs the entire nested LOSO.
#
# Trimmed to the elastic-net + plain-unpenalized path only (2026-08-24):
# a six-learner sweep (sPLS-DA, PAM, random forest, SVM) and a B-grid
# permutation-resolution sweep both existed here for the deleted
# F05/F06 screens. This repo's own prior build notes already concluded
# complex learners don't beat regularized linear models at n=16, and nested
# LOSO already needs no B-grid once B is fixed at 200 for a confirmatory
# pass rather than a screen-design comparison.

# The permutation null forks PERM_CORES workers; each would otherwise let its
# BLAS spawn one thread per core, oversubscribing the machine by orders of
# magnitude. Pin every math library to a single thread so the only parallelism
# is the fork level.
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1"
)

pacman::p_load(glmnet, pROC, parallel, withr)

ENET_ALPHAS <- c(0.25, 0.5, 0.75)
N_PERM <- 200L
PERM_SEED <- 42L
PERM_CORES <- as.integer(Sys.getenv(
  "PRED_CORES",
  as.character(max(1L, parallel::detectCores() - 2L))
))

# Train-only z-scoring: centre and scale come from the training rows and are
# applied unchanged to the held-out rows. Zero-variance columns get scale 1.
scale_train_apply <- function(x_train, x_test) {
  ctr <- colMeans(x_train)
  scl <- apply(x_train, 2, stats::sd)
  scl[scl == 0 | !is.finite(scl)] <- 1
  list(
    train = sweep(sweep(x_train, 2, ctr, "-"), 2, scl, "/"),
    test  = sweep(sweep(x_test, 2, ctr, "-"), 2, scl, "/")
  )
}

# Elastic net: inner LOSO over the training subjects tunes alpha and lambda
# together; the winning (alpha, lambda.1se) is refit on all training
# subjects. lambda.1se is the parsimonious choice; the caller can read
# lambda.min off the same fit.
fit_predict_glmnet <- function(x_tr, y_tr, x_te, family,
                               alphas = ENET_ALPHAS) {
  fold_id <- seq_len(nrow(x_tr))
  best <- NULL
  best_cvm <- Inf
  for (a in alphas) {
    cv <- cv.glmnet(x_tr, y_tr,
      family = family, alpha = a,
      foldid = fold_id, grouped = FALSE, standardize = FALSE
    )
    m <- min(cv$cvm)
    if (m < best_cvm) {
      best_cvm <- m
      best <- cv
    }
  }
  type <- if (family == "binomial") "response" else "link"
  beta <- as.matrix(coef(best, s = "lambda.1se"))[-1, 1]
  list(
    pred = as.numeric(predict(best, x_te, s = "lambda.1se", type = type)),
    selected = names(beta)[beta != 0]
  )
}

# Plain unpenalized model, only valid where p < n (the module space). Column
# names are sanitised because the trajectory space carries non-syntactic tags.
fit_predict_plain <- function(x_tr, y_tr, x_te, family) {
  df_tr <- as.data.frame(x_tr)
  df_te <- as.data.frame(x_te)
  colnames(df_tr) <- paste0("v", seq_len(ncol(df_tr)))
  colnames(df_te) <- colnames(df_tr)
  if (family == "binomial") {
    fit <- suppressWarnings(
      stats::glm(y_tr ~ ., data = df_tr, family = stats::binomial())
    )
    pred <- stats::predict(fit, df_te, type = "response")
  } else {
    fit <- stats::lm(y_tr ~ ., data = df_tr)
    pred <- stats::predict(fit, df_te)
  }
  list(pred = as.numeric(pred), selected = character(0))
}

fit_predict <- function(model, x_tr, y_tr, x_te, family) {
  switch(model,
    enet = fit_predict_glmnet(x_tr, y_tr, x_te, family),
    plain = fit_predict_plain(x_tr, y_tr, x_te, family),
    stop("unknown model: ", model)
  )
}

# One full nested LOSO pass: out-of-fold predictions aligned to y, plus the set
# of features each outer fold selected (for the selection-frequency readout).
nested_loso <- function(x, y, model, family) {
  n <- length(y)
  preds <- numeric(n)
  selected <- vector("list", n)
  for (o in seq_len(n)) {
    sc <- scale_train_apply(x[-o, , drop = FALSE], x[o, , drop = FALSE])
    fit <- fit_predict(model, sc$train, y[-o], sc$test, family)
    preds[o] <- fit$pred
    selected[[o]] <- fit$selected
  }
  list(preds = preds, selected = selected)
}

# Fraction of outer folds that selected each feature.
selection_frequency <- function(selected, model) {
  n <- length(selected)
  tab <- table(unlist(selected))
  if (!length(tab)) {
    return(data.frame(
      model = character(0), feature = character(0),
      folds = integer(0), freq = numeric(0)
    ))
  }
  data.frame(
    model = model, feature = names(tab),
    folds = as.integer(tab), freq = as.numeric(tab) / n,
    row.names = NULL
  ) |>
    dplyr::arrange(dplyr::desc(freq))
}

# Statistic from out-of-fold predictions. Class arm reports AUC; the continuous
# arm reports the leave-one-out Q^2 (higher is better).
stat_auc <- function(y, preds) {
  as.numeric(
    pROC::auc(
      pROC::roc(y, preds,
        quiet = TRUE, levels = c(0, 1),
        direction = "<"
      )
    )
  )
}

stat_q2 <- function(y, preds) {
  1 - sum((y - preds)^2) / sum((y - mean(y))^2)
}

stat_rmse <- function(y, preds) {
  sqrt(mean((y - preds)^2))
}

stat_spearman <- function(y, preds) {
  suppressWarnings(stats::cor(y, preds, method = "spearman"))
}

perm_matrix <- function(n, nperm = N_PERM, seed = PERM_SEED) {
  withr::with_seed(seed, {
    replicate(nperm, sample.int(n))
  })
}

perm_p <- function(observed, null_vals, side = c("greater", "less")) {
  side <- match.arg(side)
  null_vals <- null_vals[is.finite(null_vals)]
  hits <- if (side == "greater") {
    sum(null_vals >= observed)
  } else {
    sum(null_vals <= observed)
  }
  (hits + 1) / (length(null_vals) + 1)
}

# Run the class arm for one (feature space, model): observed AUC with a DeLong
# interval, a permutation p from re-running the whole nested LOSO, and the
# per-fold selection frequency.
run_class_cell <- function(x, y, model, nperm = N_PERM, cores = PERM_CORES) {
  fit <- nested_loso(x, y, model, "binomial")
  preds <- fit$preds
  obs <- stat_auc(y, preds)
  roc_obj <- pROC::roc(y, preds,
    quiet = TRUE, levels = c(0, 1), direction = "<"
  )
  ci <- as.numeric(pROC::ci.auc(roc_obj))

  pm <- perm_matrix(length(y), nperm)
  null_auc <- unlist(mclapply(seq_len(nperm), function(b) {
    yp <- y[pm[, b]]
    stat_auc(yp, nested_loso(x, yp, model, "binomial")$preds)
  }, mc.cores = cores))

  list(
    summary = data.frame(
      model = model, metric = "AUC", n = length(y),
      estimate = obs, ci_lo = ci[1], ci_hi = ci[3],
      perm_p = perm_p(obs, null_auc, "greater"),
      null_mean = mean(null_auc, na.rm = TRUE)
    ),
    roc = data.frame(
      model = model,
      fpr = rev(1 - roc_obj$specificities),
      tpr = rev(roc_obj$sensitivities)
    ),
    preds = data.frame(
      model = model, subject = rownames(x), y = y, pred = preds
    ),
    selection = selection_frequency(fit$selected, model)
  )
}

# Run the continuous arm for one (feature space, model, outcome): Q^2, RMSE and
# Spearman each with a permutation p, plus the per-fold selection frequency.
run_cont_cell <- function(x, y, model, outcome, nperm = N_PERM,
                          cores = PERM_CORES) {
  fit <- nested_loso(x, y, model, "gaussian")
  preds <- fit$preds
  obs_q2 <- stat_q2(y, preds)
  obs_rmse <- stat_rmse(y, preds)
  obs_rho <- stat_spearman(y, preds)

  pm <- perm_matrix(length(y), nperm)
  null <- mclapply(seq_len(nperm), function(b) {
    yp <- y[pm[, b]]
    pp <- nested_loso(x, yp, model, "gaussian")$preds
    c(
      q2 = stat_q2(yp, pp), rmse = stat_rmse(yp, pp),
      rho = stat_spearman(yp, pp)
    )
  }, mc.cores = cores)
  null <- do.call(rbind, null)

  list(
    summary = data.frame(
      outcome = outcome, model = model, n = length(y),
      q2 = obs_q2, rmse = obs_rmse, spearman = obs_rho,
      perm_p_q2 = perm_p(obs_q2, null[, "q2"], "greater"),
      perm_p_rmse = perm_p(obs_rmse, null[, "rmse"], "less"),
      perm_p_spearman = perm_p(obs_rho, null[, "rho"], "greater"),
      null_q2_mean = mean(null[, "q2"], na.rm = TRUE)
    ),
    preds = data.frame(
      outcome = outcome, model = model, subject = rownames(x),
      y = y, pred = preds
    ),
    selection = selection_frequency(fit$selected, model)
  )
}

# Align a feature matrix and an outcome vector to shared subjects with a
# non-missing outcome, and drop zero-variance features.
align_xy <- function(x, y_named) {
  subj <- intersect(rownames(x), names(y_named)[!is.na(y_named)])
  x <- x[subj, , drop = FALSE]
  y <- y_named[subj]
  vary <- apply(x, 2, function(col) stats::sd(col) > 0)
  list(x = x[, vary, drop = FALSE], y = as.numeric(y))
}
