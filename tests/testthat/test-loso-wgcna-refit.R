pacman::p_load(here, testthat)
source(here("functions", "loso_wgcna_refit.R"))

test_that("match_modules pairs by protein overlap and reports jaccard", {
  full <- c(p1 = "blue", p2 = "blue", p3 = "blue", p4 = "red", p5 = "red")
  train <- c(
    p1 = "turquoise", p2 = "turquoise", p3 = "turquoise",
    p4 = "green", p5 = "green"
  )

  m <- match_modules(train, full)

  expect_setequal(m$full, c("blue", "red"))
  expect_equal(m$train[m$full == "blue"], "turquoise")
  expect_equal(m$jaccard[m$full == "blue"], 1, tolerance = 1e-8)
  expect_equal(m$train[m$full == "red"], "green")
})

test_that("match_modules computes a partial overlap correctly", {
  full <- c(p1 = "blue", p2 = "blue", p3 = "blue", p4 = "blue")
  train <- c(p1 = "brown", p2 = "brown", p3 = "yellow", p4 = "yellow")

  m <- match_modules(train, full)

  expect_equal(m$jaccard[m$full == "blue"], 0.5, tolerance = 1e-8)
})

test_that("match_modules excludes grey from both sides", {
  full <- c(p1 = "blue", p2 = "blue", p3 = "grey")
  train <- c(p1 = "blue", p2 = "blue", p3 = "grey")

  m <- match_modules(train, full)

  expect_equal(nrow(m), 1L)
  expect_equal(m$full, "blue")
})

test_that("project_pc1 reproduces the fit on the rows it was fitted to", {
  set.seed(11)
  x <- matrix(rnorm(60), nrow = 10, dimnames = list(NULL, paste0("p", 1:6)))

  fit <- fit_pc1(x)
  direct <- project_pc1(fit, x)

  expect_length(direct, 10L)
  expect_equal(stats::sd(direct) > 0, TRUE)

  again <- project_pc1(fit_pc1(x), x)
  expect_equal(direct, again, tolerance = 1e-10)
})

test_that("project_pc1 orients the component along the average protein", {
  set.seed(12)
  base <- rnorm(10)
  x <- cbind(
    p1 = base, p2 = base + rnorm(10, sd = 0.01),
    p3 = base + rnorm(10, sd = 0.01)
  )

  score <- project_pc1(fit_pc1(x), x)

  expect_gt(stats::cor(score, rowMeans(x)), 0)
})

test_that("fit_pc1 rejects a single-protein module", {
  x <- matrix(rnorm(10), ncol = 1, dimnames = list(NULL, "p1"))
  expect_error(fit_pc1(x), "at least two")
})
