pacman::p_load(here, testthat, dplyr)
source(here("functions", "sweep_rollup.R"))

mixed_table <- function() {
  bind_rows(
    data.frame(
      level = "modules", config = "T2", outcome = "d_mcsa", model = "enet",
      B = c(0L, 200L, 1000L), q2 = 0.25, perm_p_q2 = c(NA, 0.02, 0.018)
    ),
    data.frame(
      level = "pathways", config = "T2", outcome = "d_mcsa", model = "enet",
      B = c(0L, 200L), q2 = 0.41, perm_p_q2 = c(NA, 0.015)
    ),
    data.frame(
      level = "proteins", config = "T1", outcome = "d_fcsa_I", model = "lasso",
      B = c(0L, 200L), q2 = -0.3, perm_p_q2 = c(NA, 0.8)
    )
  )
}

test_that("best_b_per_cell keeps every cell when B values differ across cells", {
  out <- best_b_per_cell(mixed_table())

  expect_equal(nrow(out), 3L)
  expect_setequal(out$level, c("modules", "pathways", "proteins"))
  expect_equal(out$B[out$level == "modules"], 1000L)
  expect_equal(out$B[out$level == "pathways"], 200L)
  expect_equal(out$B[out$level == "proteins"], 200L)
})

test_that("pooled max(B) would have dropped the lower-B cells", {
  d <- mixed_table()
  pooled <- d[d$B == max(d$B), , drop = FALSE]

  expect_equal(nrow(pooled), 1L)
  expect_equal(pooled$level, "modules")
  expect_gt(nrow(best_b_per_cell(d)), nrow(pooled))
})

test_that("cell_key adapts to a table without an outcome column", {
  classification <- data.frame(
    level = "modules", config = "total", model = "plain",
    B = c(0L, 200L), estimate = 0.78, perm_p = c(NA, 0.045)
  )

  expect_equal(cell_key(classification), c("level", "config", "model"))
  expect_equal(nrow(best_b_per_cell(classification)), 1L)
  expect_equal(best_b_per_cell(classification)$B, 200L)
})

test_that("a single-B table is returned one row per cell", {
  d <- data.frame(
    level = c("pathways", "pathways"), config = c("T1", "T2"),
    outcome = "d_mcsa", model = "spls", B = 200L, q2 = c(0.1, 0.4),
    perm_p_q2 = c(0.3, 0.01)
  )

  expect_equal(nrow(best_b_per_cell(d)), 2L)
})
