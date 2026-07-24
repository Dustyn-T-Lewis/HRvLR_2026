test_that("pred_contrast_matrix builds total and trajectory encodings", {
  source(here::here("functions", "pred_features.R"))

  feat <- matrix(
    1:18,
    nrow = 3,
    dimnames = list(
      c("g1", "g2", "g3"),
      c("s1_T1", "s1_T2", "s1_T3", "s2_T1", "s2_T2", "s2_T3")
    )
  )
  meta <- data.frame(
    sample = colnames(feat),
    subject = rep(c("s1", "s2"), each = 3),
    timepoint = factor(rep(c("T1", "T2", "T3"), 2)),
    stringsAsFactors = FALSE
  )

  total <- pred_contrast_matrix(feat, meta, "total")
  expect_equal(rownames(total), c("s1", "s2"))
  expect_equal(unname(total["s1", ]), unname(feat[, "s1_T3"] - feat[, "s1_T1"]))

  traj <- pred_contrast_matrix(feat, meta, "trajectory")
  expect_equal(ncol(traj), 9)
  expect_true(all(c("g1@T1", "g2@T2", "g3@T3") %in% colnames(traj)))
  expect_equal(traj["s2", "g1@T2"], feat["g1", "s2_T2"])
})

test_that("trajectory drops a subject missing a timepoint", {
  source(here::here("functions", "pred_features.R"))
  feat <- matrix(1:8,
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("s1_T1", "s1_T2", "s1_T3", "s2_T1"))
  )
  meta <- data.frame(
    sample = colnames(feat),
    subject = c("s1", "s1", "s1", "s2"),
    timepoint = factor(c("T1", "T2", "T3", "T1")),
    stringsAsFactors = FALSE
  )
  traj <- pred_contrast_matrix(feat, meta, "trajectory")
  expect_equal(rownames(traj), "s1")
})

test_that("sparse learners and RF recover a planted class signal", {
  source(here::here("functions", "shared_prediction.R"))
  set.seed(3)
  n <- 16L
  x <- matrix(rnorm(n * 30), n,
    dimnames = list(paste0("s", seq_len(n)), paste0("f", 1:30))
  )
  y <- rep(c(0, 1), each = 8)
  x[, "f1"] <- x[, "f1"] + y * 3

  for (m in c("enet", "lasso", "rf")) {
    res <- suppressWarnings(nested_loso(x, y, m, "binomial"))
    expect_gt(stat_auc(y, res$preds), 0.6)
  }
  # ridge and svm are breadth-only; at n = 16 under LOSO they can sit at or
  # below chance (ridge over-shrinks and inverts) as the null exposes.
  for (m in c("ridge", "svm")) {
    res <- suppressWarnings(nested_loso(x, y, m, "binomial"))
    expect_length(res$preds, n)
    expect_true(all(is.finite(res$preds)))
  }
  for (m in c("enet", "lasso")) {
    sel <- selection_frequency(
      suppressWarnings(nested_loso(x, y, m, "binomial"))$selected, m
    )
    expect_true("f1" %in% sel$feature)
  }
})

test_that("dense learners report no sparse signature", {
  source(here::here("functions", "shared_prediction.R"))
  set.seed(4)
  x <- matrix(rnorm(16 * 5), 16,
    dimnames = list(paste0("s", 1:16), paste0("f", 1:5))
  )
  for (m in c("rf", "svm")) {
    fit <- suppressWarnings(nested_loso(x, rbinom(16, 1, 0.5), m, "binomial"))
    expect_true(all(lengths(fit$selected) == 0))
  }
})

test_that("plain model fits where p < n and predicts both families", {
  source(here::here("functions", "shared_prediction.R"))
  set.seed(5)
  x <- matrix(rnorm(16 * 4), 16,
    dimnames = list(paste0("s", 1:16), paste0("f", 1:4))
  )
  yb <- as.numeric(x[, 1] > 0)
  fit_b <- suppressWarnings(nested_loso(x, yb, "plain", "binomial"))
  expect_length(fit_b$preds, 16)
  expect_true(all(fit_b$preds >= 0 & fit_b$preds <= 1))

  yc <- x[, 1] * 2 + rnorm(16, sd = 0.2)
  fit_c <- nested_loso(x, yc, "plain", "gaussian")
  expect_gt(stat_q2(yc, fit_c$preds), 0)
})

test_that("perm_matrix prefixes reproduce smaller draws", {
  source(here::here("functions", "shared_prediction.R"))
  expect_identical(perm_matrix(16, 5), perm_matrix(16, 20)[, 1:5])
})

test_that("sweep_class_cell emits one row per B and slices one null", {
  source(here::here("functions", "shared_prediction.R"))
  set.seed(6)
  x <- matrix(rnorm(16 * 8), 16,
    dimnames = list(paste0("s", 1:16), paste0("f", 1:8))
  )
  y <- rep(c(0, 1), each = 8)
  x[, "f1"] <- x[, "f1"] + y * 2

  res <- suppressWarnings(
    sweep_class_cell(x, y, "lasso", b_grid = c(0L, 10L, 20L), cores = 1L)
  )
  expect_equal(sort(res$summary$B), c(0L, 10L, 20L))
  expect_true(is.na(res$summary$perm_p[res$summary$B == 0]))
  expect_equal(nrow(res$null), 20L)
  expect_equal(length(unique(res$summary$estimate)), 1L)

  b10 <- res$summary$perm_p[res$summary$B == 10]
  hits <- sum(res$null$auc[1:10] >= res$summary$estimate[1])
  expect_equal(b10, (hits + 1) / 11)
})

test_that("sweep_cont_cell reports Q2 and a per-B permutation p", {
  source(here::here("functions", "shared_prediction.R"))
  set.seed(7)
  x <- matrix(rnorm(16 * 6), 16,
    dimnames = list(paste0("s", 1:16), paste0("f", 1:6))
  )
  y <- x[, 1] * 2 + rnorm(16, sd = 0.3)

  res <- suppressWarnings(
    sweep_cont_cell(x, y, "ridge", "d_mcsa", b_grid = c(0L, 20L), cores = 1L)
  )
  expect_equal(sort(res$summary$B), c(0L, 20L))
  expect_true(is.na(res$summary$perm_p_q2[res$summary$B == 0]))
  expect_equal(nrow(res$null), 20L)
  expect_true(all(res$summary$outcome == "d_mcsa"))
})
