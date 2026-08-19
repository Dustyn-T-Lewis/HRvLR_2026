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

test_that("the size floor counts detected members, not annotated ones", {
  source(here::here("functions", "shared_singscore.R"))

  set.seed(2)
  expr <- matrix(rnorm(50 * 6),
    nrow = 50,
    dimnames = list(paste0("g", 1:50), paste0("s", 1:6))
  )
  # big_but_absent is annotated 200-wide and measured 3-wide: the case that
  # motivated the floor. measured is 20 of 20.
  sets <- list(
    big_but_absent = c(paste0("g", 1:3), paste0("absent", 1:197)),
    measured = paste0("g", 21:40)
  )

  kept <- score_singscore(expr, sets, min_size = 15L)
  expect_identical(rownames(kept), "measured")
  expect_equal(nrow(kept), 1L)

  expect_setequal(
    rownames(suppressWarnings(score_singscore(expr, sets, min_size = 1L))),
    c("big_but_absent", "measured")
  )
})
