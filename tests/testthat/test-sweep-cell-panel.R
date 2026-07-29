pacman::p_load(here, testthat)
source(here("functions", "sweep_grid.R"))
source(here("functions", "sweep_cell_panel.R"))

test_that("SCREEN_SIZE totals match the grid definitions", {
  assoc_per_leaf <- length(ASSOC_CONT_METHODS) * length(ADAPT_OUTCOMES) +
    length(ASSOC_GROUP_METHODS)
  assoc_total <- length(SWEEP_LEVELS) * length(SWEEP_CONFIGS) * assoc_per_leaf

  grid_rows <- function(methods_grid) {
    sum(vapply(SWEEP_LEVELS, function(lv) {
      methods <- methods_for_level(methods_grid, lv)
      n <- length(methods) * length(SWEEP_CONFIGS)
      n - as.integer("plain" %in% methods)
    }, integer(1)))
  }

  expect_equal(assoc_total, SCREEN_SIZE[["F04_association"]])
  expect_equal(grid_rows(CLASS_METHODS), SCREEN_SIZE[["F05_classification"]])
  expect_equal(
    grid_rows(CONT_METHODS) * length(ADAPT_OUTCOMES),
    SCREEN_SIZE[["F06_prediction"]]
  )

  expect_equal(unname(SCREEN_SIZE[["F04_association"]]), 420L)
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

test_that("wilcoxon-shaped input scores on effect, not an all-NA t", {
  cell <- data.frame(
    feature = paste0("PROT", seq_len(5)),
    t = NA_real_,
    effect = c(0.4, -0.9, 0.1, -0.3, 0.6),
    p = c(0.01, 0.2, 0.5, 0.3, 0.04)
  )

  sc <- assoc_score(cell)

  expect_identical(sc$axis, "effect")
  expect_identical(sc$stat, "top absolute effect")
  expect_equal(sc$value, cell$effect)

  footer_stat <- max(abs(sc$value), na.rm = TRUE)
  expect_true(is.finite(footer_stat))
})
