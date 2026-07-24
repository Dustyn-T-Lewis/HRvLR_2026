# Leakage-free nested leave-one-subject-out harness for the prediction suite.
# Rows are subjects, so leave-one-subject-out is leave-one-row-out. Every fold:
# features are z-scored on the training subjects and that centre/scale is
# applied to the held-out subject; each engine's hyperparameters are tuned by an
# inner LOSO on the training subjects only; the held-out subject is never seen
# until it is predicted. Every fit also reports which features it selected, so
# the outer folds yield a per-feature selection frequency (methods doc Part 6).
# The permutation null shuffles the outcome across subjects and re-runs the
# entire nested LOSO.

pacman::p_load(glmnet, mixOmics, pamr, pROC, ranger, e1071, parallel, withr)

GLMNET_ALPHAS <- c(0.1, 0.5, 1.0)
ENET_ALPHAS <- c(0.25, 0.5, 0.75)
RF_TREES <- 500L
SVM_COST <- 1
SPLS_NCOMP <- 1L
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

# keepX candidates scale with the feature count: a handful for the low-dimension
# module space, a wider sparse sweep for pathways and proteins. Capped at p so a
# 10-feature space never asks for 50.
spls_keepx_grid <- function(p) {
  unique(pmin(c(2L, 5L, 10L, 20L, 50L), p))
}

# glmnet: inner LOSO over the training subjects tunes alpha and lambda together;
# the winning (alpha, lambda.1se) is refit on all training subjects. lambda.1se
# is the parsimonious choice (methods doc Part 4); the caller can read
# lambda.min off the same fit.
fit_predict_glmnet <- function(x_tr, y_tr, x_te, family,
                               alphas = GLMNET_ALPHAS) {
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

spls_score <- function(fit, x_te, family) {
  pr <- predict(fit, x_te)
  if (family == "binomial") {
    pr$predict[, "1", SPLS_NCOMP]
  } else {
    pr$predict[, 1, SPLS_NCOMP]
  }
}

# sPLS / sPLS-DA: inner LOSO over the training subjects tunes keepX; the winning
# keepX is refit on all training subjects. Scaling is off because the harness
# already z-scored on train.
fit_predict_spls <- function(x_tr, y_tr, x_te, family) {
  grid <- spls_keepx_grid(ncol(x_tr))
  fit_one <- function(xx, yy, keepx) {
    if (family == "binomial") {
      splsda(xx, factor(yy, levels = c(0, 1)),
        ncomp = SPLS_NCOMP, keepX = keepx, scale = FALSE
      )
    } else {
      spls(xx, yy, ncomp = SPLS_NCOMP, keepX = keepx, scale = FALSE)
    }
  }
  best_k <- grid[[1]]
  best_score <- Inf
  for (k in grid) {
    err <- 0
    for (i in seq_len(nrow(x_tr))) {
      sc <- scale_train_apply(x_tr[-i, , drop = FALSE], x_tr[i, , drop = FALSE])
      fit <- fit_one(sc$train, y_tr[-i], k)
      pred <- spls_score(fit, sc$test, family)
      err <- err + (y_tr[i] - pred)^2
    }
    if (err < best_score) {
      best_score <- err
      best_k <- k
    }
  }
  fit <- fit_one(x_tr, y_tr, best_k)
  sel <- mixOmics::selectVar(fit, comp = SPLS_NCOMP)
  list(
    pred = spls_score(fit, x_te, family),
    selected = if (family == "binomial") sel$name else sel$X$name
  )
}

# Nearest shrunken centroids (Tibshirani 2002): the interpretable diagonal
# baseline the small-n literature favours. pamr wants features x samples, so the
# harness matrices are transposed on the way in. The shrinkage threshold is
# tuned by a leave-one-out pamr.cv on training subjects; ties break to largest
# threshold (sparsest centroid). Classification only.
fit_predict_pam <- function(x_tr, y_tr, x_te) {
  dat <- list(
    x = t(x_tr), y = factor(y_tr, levels = c(0, 1)),
    geneid = colnames(x_tr)
  )
  utils::capture.output(fit <- pamr::pamr.train(dat))
  folds <- as.list(seq_len(nrow(x_tr)))
  utils::capture.output(cv <- pamr::pamr.cv(fit, dat, folds = folds))
  thr <- max(cv$threshold[cv$error == min(cv$error)])
  post <- pamr::pamr.predict(fit, t(x_te), threshold = thr, type = "posterior")
  nz <- pamr::pamr.predict(fit, t(x_te), threshold = thr, type = "nonzero")
  list(pred = post[, "1"], selected = colnames(x_tr)[nz])
}

# Random forest (ranger): a breadth-only robustness learner, single-threaded so
# it nests cleanly inside the permutation mclapply, seeded for reproducibility.
# It uses every feature, so it yields no sparse signature; selected is empty by
# design (the methods doc forbids ranking off RF importance).
fit_predict_ranger <- function(x_tr, y_tr, x_te, family) {
  df_tr <- as.data.frame(x_tr)
  df_te <- as.data.frame(x_te)
  if (family == "binomial") {
    rf <- ranger::ranger(
      x = df_tr, y = factor(y_tr, levels = c(0, 1)),
      probability = TRUE, num.trees = RF_TREES, num.threads = 1L, seed = 1L
    )
    pred <- predict(rf, df_te)$predictions[, "1"]
  } else {
    rf <- ranger::ranger(
      x = df_tr, y = y_tr,
      num.trees = RF_TREES, num.threads = 1L, seed = 1L
    )
    pred <- predict(rf, df_te)$predictions
  }
  list(pred = as.numeric(pred), selected = character(0))
}

# Linear-kernel SVM (e1071) at a fixed cost: nothing is tuned, so there is no
# in-fold choice to leak. A dense robustness learner, no sparse signature.
fit_predict_svm <- function(x_tr, y_tr, x_te, family) {
  if (family == "binomial") {
    fit <- e1071::svm(x_tr, factor(y_tr, levels = c(0, 1)),
      kernel = "linear", cost = SVM_COST, scale = FALSE, probability = TRUE
    )
    pr <- predict(fit, x_te, probability = TRUE)
    pred <- attr(pr, "probabilities")[, "1"]
  } else {
    fit <- e1071::svm(x_tr, y_tr,
      kernel = "linear", cost = SVM_COST, scale = FALSE
    )
    pred <- predict(fit, x_te)
  }
  list(pred = as.numeric(pred), selected = character(0))
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
    glmnet = fit_predict_glmnet(x_tr, y_tr, x_te, family),
    enet = fit_predict_glmnet(x_tr, y_tr, x_te, family, ENET_ALPHAS),
    lasso = fit_predict_glmnet(x_tr, y_tr, x_te, family, 1),
    ridge = fit_predict_glmnet(x_tr, y_tr, x_te, family, 0),
    spls = fit_predict_spls(x_tr, y_tr, x_te, family),
    pam = fit_predict_pam(x_tr, y_tr, x_te),
    rf = fit_predict_ranger(x_tr, y_tr, x_te, family),
    svm = fit_predict_svm(x_tr, y_tr, x_te, family),
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
  as.numeric(pROC::auc(pROC::roc(y, preds,
    quiet = TRUE, levels = c(0, 1),
    direction = "<"
  )))
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

B_GRID <- c(0L, 200L, 1000L)

# Metric null under B outcome shuffles, each re-running the whole nested LOSO.
# Because perm_matrix is seeded, the first b columns of a bmax draw reproduce a
# standalone b-permutation null, so one bmax pass serves every smaller B in the
# grid by prefix-slicing.
sweep_class_null <- function(x, y, model, bmax, cores) {
  pm <- perm_matrix(length(y), bmax)
  unlist(mclapply(seq_len(bmax), function(b) {
    yp <- y[pm[, b]]
    stat_auc(yp, nested_loso(x, yp, model, "binomial")$preds)
  }, mc.cores = cores))
}

sweep_cont_null <- function(x, y, model, bmax, cores) {
  pm <- perm_matrix(length(y), bmax)
  null <- mclapply(seq_len(bmax), function(b) {
    yp <- y[pm[, b]]
    pp <- nested_loso(x, yp, model, "gaussian")$preds
    c(q2 = stat_q2(yp, pp), rho = stat_spearman(yp, pp))
  }, mc.cores = cores)
  do.call(rbind, null)
}

# Classification cell across a grid of permutation counts: observed AUC once,
# one null pass at max(B), one summary row per B with that B's permutation p.
sweep_class_cell <- function(x, y, model, b_grid = B_GRID, cores = PERM_CORES) {
  fit <- nested_loso(x, y, model, "binomial")
  preds <- fit$preds
  obs <- stat_auc(y, preds)
  roc_obj <- pROC::roc(y, preds,
    quiet = TRUE, levels = c(0, 1), direction = "<"
  )
  ci <- as.numeric(pROC::ci.auc(roc_obj))
  bmax <- max(b_grid)
  null <- if (bmax > 0) {
    sweep_class_null(x, y, model, bmax, cores)
  } else {
    numeric(0)
  }
  summ <- do.call(rbind, lapply(sort(b_grid), function(b) {
    nb <- if (b > 0) null[seq_len(b)] else numeric(0)
    data.frame(
      model = model, metric = "AUC", n = length(y), p = ncol(x), B = b,
      estimate = obs, ci_lo = ci[1], ci_hi = ci[3],
      perm_p = if (b > 0) perm_p(obs, nb, "greater") else NA_real_,
      null_mean = if (b > 0) mean(nb) else NA_real_,
      null_sd = if (b > 0) stats::sd(nb) else NA_real_
    )
  }))
  list(
    summary = summ, null = data.frame(model = model, auc = null),
    preds = data.frame(
      model = model, subject = rownames(x), y = y, pred = preds
    ),
    selection = selection_frequency(fit$selected, model)
  )
}

# Continuous cell for one outcome across a grid of permutation counts.
sweep_cont_cell <- function(x, y, model, outcome, b_grid = B_GRID,
                            cores = PERM_CORES) {
  fit <- nested_loso(x, y, model, "gaussian")
  preds <- fit$preds
  obs_q2 <- stat_q2(y, preds)
  obs_rho <- stat_spearman(y, preds)
  bmax <- max(b_grid)
  null <- if (bmax > 0) sweep_cont_null(x, y, model, bmax, cores) else NULL
  summ <- do.call(rbind, lapply(sort(b_grid), function(b) {
    nq <- if (b > 0) null[seq_len(b), "q2"] else numeric(0)
    nr <- if (b > 0) null[seq_len(b), "rho"] else numeric(0)
    data.frame(
      outcome = outcome, model = model, n = length(y), p = ncol(x), B = b,
      q2 = obs_q2, spearman = obs_rho,
      perm_p_q2 = if (b > 0) perm_p(obs_q2, nq, "greater") else NA_real_,
      perm_p_spearman = if (b > 0) perm_p(obs_rho, nr, "greater") else NA_real_,
      null_q2_mean = if (b > 0) mean(nq) else NA_real_,
      null_q2_sd = if (b > 0) stats::sd(nq) else NA_real_
    )
  }))
  null_df <- if (is.null(null)) {
    data.frame(outcome = outcome, model = model, q2 = numeric(0))
  } else {
    data.frame(outcome = outcome, model = model, q2 = null[, "q2"])
  }
  list(
    summary = summ, null = null_df,
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
