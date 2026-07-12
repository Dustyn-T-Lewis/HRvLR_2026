pacman::p_load(here, testthat)
source(here("05_test", "modules", "a_script", "functions", "honest_refit.R"))

test_that("jaccard and pc1 projection behave", {
  expect_equal(jaccard(c("a", "b", "c"), c("b", "c", "d")), 2 / 4)

  set.seed(3)
  m <- matrix(rnorm(30 * 6), nrow = 30)
  fit <- fit_pc1(m[1:24, ])
  s_in <- project_pc1(fit, m[1:24, ])
  s_out <- project_pc1(fit, m[25:30, ])
  expect_length(s_out, 6)
  expect_gt(abs(cor(s_in, rowMeans(scale(m[1:24, ])))), 0.5)
})

test_that("match_modules pairs by best jaccard", {
  full <- c(p1 = "red", p2 = "red", p3 = "blue", p4 = "blue")
  train <- c(p1 = "green", p2 = "green", p3 = "yellow", p4 = "yellow")
  mm <- match_modules(full, train)
  expect_equal(mm$train[mm$full == "red"], "green")
  expect_equal(mm$jaccard[mm$full == "red"], 1)
})
