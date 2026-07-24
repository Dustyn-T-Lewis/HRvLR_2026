test_that("shared_prediction harness is deterministic and leakage-free", {
  source(here::here("functions", "shared_prediction.R"))

  set.seed(1)
  x <- matrix(rnorm(20 * 6),
    nrow = 20, dimnames = list(paste0("s", 1:20), NULL)
  )
  y <- as.numeric(x[, 1] > 0)

  fit1 <- suppressWarnings(nested_loso(x, y, "glmnet", "binomial"))
  fit2 <- suppressWarnings(nested_loso(x, y, "glmnet", "binomial"))
  expect_identical(fit1$preds, fit2$preds)
  expect_length(fit1$preds, 20)
  expect_length(fit1$selected, 20)

  # perm_matrix must be reproducible from its seed alone
  expect_identical(perm_matrix(20, nperm = 5), perm_matrix(20, nperm = 5))

  # align_xy drops zero-variance columns and NA outcomes
  xz <- cbind(x, const = 1)
  al <- align_xy(xz, setNames(c(y[1:19], NA), rownames(xz)))
  expect_false("const" %in% colnames(al$x))
  expect_equal(nrow(al$x), 19)
})

test_that("all three model paths recover a planted signal and report it", {
  source(here::here("functions", "shared_prediction.R"))

  set.seed(2)
  n <- 16L
  x <- matrix(rnorm(n * 30), n,
    dimnames = list(paste0("s", seq_len(n)), paste0("f", 1:30))
  )
  y <- rep(c(0, 1), each = 8)
  x[, "f1"] <- x[, "f1"] + y * 3

  for (m in c("glmnet", "spls", "pam")) {
    res <- suppressWarnings(run_class_cell(x, y, m, nperm = 20, cores = 1))
    expect_gt(res$summary$estimate, 0.5)
    expect_true("f1" %in% res$selection$feature)
    expect_true(all(res$selection$freq >= 0 & res$selection$freq <= 1))
  }
})

test_that("selection_frequency counts folds and stays in [0, 1]", {
  source(here::here("functions", "shared_prediction.R"))

  sel <- list(c("a", "b"), c("a"), c("a", "c"), character(0))
  freq <- selection_frequency(sel, "glmnet")
  expect_equal(freq$freq[freq$feature == "a"], 0.75)
  expect_equal(freq$folds[freq$feature == "b"], 1L)
  expect_true(all(freq$freq <= 1))

  empty <- selection_frequency(list(character(0), character(0)), "pam")
  expect_equal(nrow(empty), 0L)
})

test_that("spls keepX grid is data-driven and capped at p", {
  source(here::here("functions", "shared_prediction.R"))
  expect_true(all(spls_keepx_grid(6) <= 6))
  expect_setequal(spls_keepx_grid(6), c(2, 5, 6))
  expect_true(max(spls_keepx_grid(2000)) == 50)
})

test_that("stat_q2 is 1 for perfect fit and <=0 for mean prediction", {
  source(here::here("functions", "shared_prediction.R"))
  y <- c(1, 2, 3, 4, 5)
  expect_equal(stat_q2(y, y), 1, tolerance = 1e-12)
  expect_equal(stat_q2(y, rep(mean(y), 5)), 0, tolerance = 1e-12)
  expect_lt(stat_q2(y, rev(y)), 0)
})

test_that("perm_p respects the 1/(B+1) floor and side", {
  source(here::here("functions", "shared_prediction.R"))
  expect_equal(perm_p(1, rep(0, 9), "greater"), 1 / 10)
  expect_equal(perm_p(0, rep(1, 9), "less"), 1 / 10)
  expect_equal(perm_p(0.5, c(1, 1, 0, 0), "greater"), 3 / 5)
})
