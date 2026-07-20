test_that("score_singscore matches per-sample simpleScore and is leakage-free", {
  source(here::here("functions", "shared_singscore.R"))
  pacman::p_load(singscore)

  set.seed(1)
  expr <- matrix(rnorm(50 * 8),
    nrow = 50,
    dimnames = list(paste0("g", 1:50), paste0("s", 1:8))
  )
  sets <- list(A = paste0("g", 1:10), B = paste0("g", 11:25))

  got <- score_singscore(expr, sets, min_size = 1L)
  ranks <- rankGenes(expr)
  want_A <- simpleScore(ranks, upSet = sets$A)$TotalScore
  expect_equal(unname(got["A", ]), want_A, tolerance = 1e-10)
  expect_equal(dim(got), c(2, 8))

  # leakage-free: scoring a subset of columns gives the same per-sample score
  sub <- score_singscore(expr[, 1:4], sets, min_size = 1L)
  expect_equal(got[, 1:4], sub, tolerance = 1e-10)
})
