pacman::p_load(here, testthat)
source(here("functions", "sweep_grid.R"))
source(here("functions", "sweep_cell_panel.R"))

test_that("SCREEN_SIZE totals match the grid definitions", {
  grid_rows <- function(methods_grid) {
    sum(vapply(SWEEP_LEVELS, function(lv) {
      methods <- methods_for_level(methods_grid, lv)
      n <- length(methods) * length(SWEEP_CONFIGS)
      n - as.integer("plain" %in% methods)
    }, integer(1)))
  }

  expect_equal(grid_rows(CLASS_METHODS), SCREEN_SIZE[["F05_classification"]])
  expect_equal(
    grid_rows(CONT_METHODS) * length(ADAPT_OUTCOMES),
    SCREEN_SIZE[["F06_prediction"]]
  )

  expect_equal(unname(SCREEN_SIZE[["F05_classification"]]), 153L)
  expect_equal(unname(SCREEN_SIZE[["F06_prediction"]]), 792L)
})

test_that("ridge shows no drivers while lasso does", {
  selection <- data.frame(
    feature = paste0("PROT", seq_len(12)),
    score = seq(0.1, 1.2, length.out = 12)
  )

  p_ridge <- drivers_or_empty(
    selection, "proteins", "ridge", "fold selection freq.", NULL, "drivers"
  )
  p_lasso <- drivers_or_empty(
    selection, "proteins", "lasso", "fold selection freq.", NULL, "drivers"
  )

  expect_true(is.null(nrow(p_ridge$data)))
  expect_false(is.null(nrow(p_lasso$data)))
  expect_true(nrow(p_lasso$data) > 0)
})
