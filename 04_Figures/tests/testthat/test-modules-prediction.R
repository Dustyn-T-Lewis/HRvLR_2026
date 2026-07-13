pacman::p_load(here, testthat)
source(here("04_Figures", "05_test", "modules", "a_script", "functions", "prediction.R"))

test_that("module_auc is 1 for separable data and perm p is small", {
  y <- c(rep(0, 8), rep(1, 8))
  x <- c(rnorm(8, 0), rnorm(8, 6))
  expect_equal(module_auc(y, x), 1, tolerance = 1e-8)
  expect_lt(perm_p_unpaired(y, x, n_perm = 200), 0.05)
})

test_that("perm p is near uniform for noise", {
  set.seed(42)
  y <- rep(c(0, 1), each = 8)
  x <- rnorm(16)
  expect_gt(perm_p_unpaired(y, x, n_perm = 200), 0.05)
})
