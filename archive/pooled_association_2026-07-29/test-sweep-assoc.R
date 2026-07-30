test_that("assoc_lm matches base lm slope and p per feature", {
  source(here::here("functions", "sweep_assoc.R"))
  set.seed(11)
  x <- matrix(rnorm(16 * 3), 16, dimnames = list(NULL, c("a", "b", "c")))
  y <- 2 * x[, "a"] + rnorm(16, sd = 0.5)

  out <- assoc_lm(x, y)
  ref <- summary(lm(y ~ x[, "a"]))$coefficients[2, ]
  expect_equal(out$effect[out$feature == "a"], unname(ref["Estimate"]),
    tolerance = 1e-8
  )
  expect_equal(out$p[out$feature == "a"], unname(ref["Pr(>|t|)"]),
    tolerance = 1e-8
  )
  expect_true(all(out$bh >= 0 & out$bh <= 1))
})

test_that("assoc_spearman matches cor.test rho and approx p", {
  source(here::here("functions", "sweep_assoc.R"))
  set.seed(12)
  x <- matrix(rnorm(16 * 2), 16, dimnames = list(NULL, c("a", "b")))
  y <- rank(x[, "a"]) + rnorm(16)

  out <- assoc_spearman(x, y)
  ct <- suppressWarnings(cor.test(x[, "a"], y, method = "spearman"))
  expect_equal(out$effect[out$feature == "a"], unname(ct$estimate),
    tolerance = 1e-8
  )
  expect_equal(out$p[out$feature == "a"], ct$p.value, tolerance = 0.02)
})

test_that("assoc_limma_cont recovers the slope sign and shape", {
  source(here::here("functions", "sweep_assoc.R"))
  set.seed(13)
  x <- matrix(rnorm(16 * 5), 16, dimnames = list(NULL, paste0("f", 1:5)))
  y <- 3 * x[, 1] + rnorm(16, sd = 0.3)
  out <- assoc_limma_cont(x, y)
  expect_named(out, c("feature", "effect", "t", "p", "bh"))
  expect_gt(out$effect[out$feature == "f1"], 0)
  expect_lt(out$p[out$feature == "f1"], 0.05)
})

test_that("group methods orient the effect as HR minus LR", {
  source(here::here("functions", "sweep_assoc.R"))
  set.seed(14)
  x <- matrix(rnorm(16 * 4), 16, dimnames = list(NULL, paste0("f", 1:4)))
  y <- rep(c(0, 1), each = 8)
  x[, 1] <- x[, 1] + y * 4

  lim <- assoc_limma_group(x, y)
  wil <- assoc_wilcox_group(x, y)
  expect_gt(lim$effect[lim$feature == "f1"], 0)
  expect_gt(wil$effect[wil$feature == "f1"], 0)
  expect_lt(wil$p[wil$feature == "f1"], 0.05)
  expect_true(all(is.na(wil$t)))
})

test_that("run_assoc_cell dispatches and rejects unknown methods", {
  source(here::here("functions", "sweep_assoc.R"))
  set.seed(15)
  x <- matrix(rnorm(16 * 3), 16, dimnames = list(NULL, c("a", "b", "c")))
  y_cont <- x[, "a"] + rnorm(16)
  y_grp <- rep(c(0, 1), each = 8)

  expect_equal(
    run_assoc_cell(x, y_cont, "lm", FALSE)$effect,
    assoc_lm(x, y_cont)$effect
  )
  expect_s3_class(run_assoc_cell(x, y_grp, "wilcoxon", TRUE), "tbl_df")
  expect_error(run_assoc_cell(x, y_cont, "wilcoxon", FALSE))
  expect_error(run_assoc_cell(x, y_grp, "lm", TRUE))
})
