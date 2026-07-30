pacman::p_load(here, testthat, withr, openxlsx)
source(here("functions", "sweep_grid.R"))

fake_bundle <- function(n_prot = 20, seed = 5) {
  set.seed(seed)
  m <- matrix(rnorm(n_prot * 6), nrow = n_prot)
  rownames(m) <- paste0("p", seq_len(n_prot))
  colnames(m) <- paste0("s", 1:6)
  list(feature_sets = list(proteins = m, singscore = m[1:5, ]))
}

test_that("the fingerprint is stable for identical input", {
  expect_identical(
    sweep_fingerprint(fake_bundle()), sweep_fingerprint(fake_bundle())
  )
})

test_that("the fingerprint moves when the feature set changes", {
  a <- sweep_fingerprint(fake_bundle(n_prot = 20))
  b <- sweep_fingerprint(fake_bundle(n_prot = 19))
  expect_false(identical(a, b))
})

test_that("the fingerprint moves when values change but shape does not", {
  a <- sweep_fingerprint(fake_bundle(seed = 1))
  b <- sweep_fingerprint(fake_bundle(seed = 2))
  expect_false(identical(a, b))
})

test_that("a leaf with no workbook is not done", {
  root <- withr::local_tempdir()
  expect_false(
    leaf_done("X", "proteins", "total", "rf", "abc", root_dir = root)
  )
})

test_that("a leaf is done only when its fingerprint matches", {
  root <- withr::local_tempdir()
  d <- file.path(root, "proteins", "total", "rf", "c_data")
  dir.create(d, recursive = TRUE)
  write_sweep_workbook(
    file.path(d, "results.xlsx"), list(metrics = data.frame(q2 = 1)),
    fingerprint = "abc"
  )

  expect_true(
    leaf_done("X", "proteins", "total", "rf", "abc", root_dir = root)
  )
  expect_false(
    leaf_done("X", "proteins", "total", "rf", "zzz", root_dir = root)
  )
})

test_that("a workbook written without a fingerprint is never treated as done", {
  root <- withr::local_tempdir()
  d <- file.path(root, "proteins", "total", "rf", "c_data")
  dir.create(d, recursive = TRUE)
  write.xlsx(list(metrics = data.frame(q2 = 1)), file.path(d, "results.xlsx"))

  expect_false(
    leaf_done("X", "proteins", "total", "rf", "abc", root_dir = root)
  )
})

test_that("write_sweep_workbook keeps the sheets it was given", {
  root <- withr::local_tempdir()
  p <- file.path(root, "results.xlsx")
  write_sweep_workbook(
    p, list(metrics = data.frame(q2 = 1), folds = data.frame(i = 1:3)),
    fingerprint = "abc"
  )
  expect_true(all(c("metrics", "folds") %in% getSheetNames(p)))
  expect_equal(read.xlsx(p, "folds")$i, 1:3)
})
