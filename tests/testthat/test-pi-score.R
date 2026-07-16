pacman::p_load(testthat, here, dplyr, tibble, readr)
source(here::here("03_DEP", "contrasts.R"))

test_that("pi_score is p raised to the absolute log2 fold change", {
  res <- tibble(P.Value = c(0.05, 0.5, 0.001), logFC = c(2, -3, 0.5))
  out <- add_pi_score(res)
  expect_equal(out$pi_score, c(0.05^2, 0.5^3, 0.001^0.5))
})

test_that("sig_pi carries the direction of a gated protein", {
  res <- tibble(
    P.Value = c(0.001, 0.001, 0.9, 0.9),
    logFC = c(3, -3, 3, -3)
  )
  expect_equal(add_pi_score(res)$sig_pi, c(1L, -1L, 0L, 0L))
})

test_that("the gate is strict and reads the threshold argument", {
  # pi_score is exactly 0.05 at p = 0.05, logFC = 1; `<` must exclude it.
  res <- tibble(P.Value = 0.05, logFC = 1)
  expect_equal(add_pi_score(res, pi_thresh = 0.05)$sig_pi, 0L)
  expect_equal(add_pi_score(res, pi_thresh = 0.06)$sig_pi, 1L)
})

test_that("a zero fold change is never selected", {
  # p^0 == 1 for any p, so the gate can never open on logFC == 0.
  res <- tibble(P.Value = c(1e-12, 0.5), logFC = c(0, 0))
  out <- add_pi_score(res)
  expect_equal(out$pi_score, c(1, 1))
  expect_equal(out$sig_pi, c(0L, 0L))
})

test_that("untested proteins land on 0, not NA", {
  # limma returns NA for the 34 proteins with an inestimable cell mean. The final
  # case_when branch maps them to 0L; sig_pi must never be NA.
  res <- tibble(P.Value = c(NA_real_, 0.001), logFC = c(NA_real_, 2))
  out <- add_pi_score(res)
  expect_true(is.na(out$pi_score[1]))
  expect_equal(out$sig_pi, c(0L, 1L))
  expect_false(anyNA(out$sig_pi))
})

test_that("add_pi_score reproduces the shipped DEP table exactly", {
  # Equivalence check against a table written by the inline code this function
  # replaced: recomputing pi from its own P.Value and logFC must return the
  # pi_score and sig_pi it already carries, or the refactor moved a published number.
  csv <- here::here(
    "03_DEP", "b_imputed", "c_data", "missforest", "combined_results_pi.csv"
  )
  skip_if_not(file.exists(csv), "run 03_DEP/b_imputed first")
  shipped <- read_csv(csv, show_col_types = FALSE)

  recomputed <- add_pi_score(select(shipped, P.Value, logFC))
  expect_equal(recomputed$pi_score, shipped$pi_score)
  expect_equal(recomputed$sig_pi, as.integer(shipped$sig_pi))
})
