# Leakage-free nested leave-one-subject-out harness for the prediction suite.
# Rows are subjects, so leave-one-subject-out is leave-one-row-out. Every fold:
# features are z-scored on the training subjects and that centre/scale is
# applied to the held-out subject; glmnet (alpha/lambda) and sPLS (keepX) are
# tuned by an inner LOSO on the training subjects only; the held-out subject is
# never seen until it is predicted. The permutation null shuffles the outcome
# across subjects and re-runs the entire nested LOSO.

pacman::p_load(glmnet, mixOmics, pROC, parallel)

GLMNET_ALPHAS <- c(0.1, 0.5, 1.0)
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

spls_keepx_grid <- function(p) {
  unique(pmin(c(10L, 30L), p))
}

# glmnet: inner LOSO over the training subjects tunes alpha and lambda together;
# the winning (alpha, lambda.min) is refit on all training subjects.
fit_predict_glmnet <- function(x_tr, y_tr, x_te, family) {
  fold_id <- seq_len(nrow(x_tr))
  best <- NULL
  best_cvm <- Inf
  for (a in GLMNET_ALPHAS) {
    cv <- cv.glmnet(x_tr, y_tr,
      family = family, alpha = a,
      foldid = fold_id, grouped = FALSE, standardize = FALSE
    )
    m <- min(cv$cvm)
    if (m < best_cvm) {
      best_cvm <- m
      best <- list(cv = cv)
    }
  }
  type <- if (family == "binomial") "response" else "link"
  as.numeric(predict(best$cv, x_te, s = "lambda.min", type = type))
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
  spls_score(fit, x_te, family)
}

# One full nested LOSO pass: returns out-of-fold predictions aligned to y.
nested_loso <- function(x, y, model, family) {
  n <- length(y)
  preds <- numeric(n)
  for (o in seq_len(n)) {
    sc <- scale_train_apply(x[-o, , drop = FALSE], x[o, , drop = FALSE])
    preds[o] <- if (model == "glmnet") {
      fit_predict_glmnet(sc$train, y[-o], sc$test, family)
    } else {
      fit_predict_spls(sc$train, y[-o], sc$test, family)
    }
  }
  preds
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
  hits <- if (side == "greater") sum(null_vals >= observed) else sum(null_vals <= observed)
  (hits + 1) / (length(null_vals) + 1)
}

# Run the class arm for one (feature space, model): observed AUC with a DeLong
# CV interval and a permutation p from re-running the whole nested LOSO.
run_class_cell <- function(x, y, model, nperm = N_PERM, cores = PERM_CORES) {
  preds <- nested_loso(x, y, model, "binomial")
  obs <- stat_auc(y, preds)
  ci <- as.numeric(pROC::ci.auc(pROC::roc(y, preds,
    quiet = TRUE,
    levels = c(0, 1), direction = "<"
  )))
  roc_obj <- pROC::roc(y, preds, quiet = TRUE, levels = c(0, 1), direction = "<")

  pm <- perm_matrix(length(y), nperm)
  null_auc <- unlist(mclapply(seq_len(nperm), function(b) {
    yp <- y[pm[, b]]
    stat_auc(yp, nested_loso(x, yp, model, "binomial"))
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
    preds = data.frame(model = model, subject = rownames(x), y = y, pred = preds)
  )
}

# Run the continuous arm for one (feature space, model, outcome): Q^2, RMSE and
# Spearman, each with a permutation p.
run_cont_cell <- function(x, y, model, outcome, nperm = N_PERM,
                          cores = PERM_CORES) {
  preds <- nested_loso(x, y, model, "gaussian")
  obs_q2 <- stat_q2(y, preds)
  obs_rmse <- stat_rmse(y, preds)
  obs_rho <- stat_spearman(y, preds)

  pm <- perm_matrix(length(y), nperm)
  null <- mclapply(seq_len(nperm), function(b) {
    yp <- y[pm[, b]]
    pp <- nested_loso(x, yp, model, "gaussian")
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
    )
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
