pacman::p_load(here, testthat, dplyr)
source(here("functions", "pred_complete_case.R"))

test_that("the complete-case matrix contains no missing value", {
  cc <- complete_case_matrix()
  expect_false(anyNA(cc$mat))
  expect_identical(length(cc$gene), nrow(cc$mat))
})

test_that("it is a strict subset of the normalised matrix, not a re-read", {
  cc <- complete_case_matrix()
  full <- read_csv(here("02_Normalization", "c_data", "normalized.csv"),
    show_col_types = FALSE
  )
  expect_lt(nrow(cc$mat), nrow(full))
  expect_true(all(rownames(cc$mat) %in% full$uniprot_id))
  expect_identical(ncol(cc$mat), ncol(full) - 4L)
})

test_that("the bundle carries no imputed values and no eigengenes", {
  b <- pred_load_complete_case()

  expect_setequal(names(b$feature_sets), c("proteins", "singscore"))
  expect_false(anyNA(b$feature_sets$proteins))
  expect_false(anyNA(b$feature_sets$singscore))
  expect_gt(nrow(b$feature_sets$proteins), 0L)
  expect_gt(nrow(b$feature_sets$singscore), 0L)
})

test_that("its samples and subjects match the imputed bundle exactly", {
  b <- pred_load_complete_case()
  ref <- pred_load()

  expect_setequal(colnames(b$feature_sets$proteins), ref$meta$sample)
  expect_identical(nrow(b$meta), nrow(ref$meta))
  expect_setequal(b$meta$subject, ref$meta$subject)
})

test_that("it holds strictly fewer proteins than the imputed bundle", {
  b <- pred_load_complete_case()
  ref <- pred_load()
  expect_lt(
    nrow(b$feature_sets$proteins), nrow(ref$feature_sets$proteins)
  )
})
