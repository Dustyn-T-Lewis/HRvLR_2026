fixture_meta <- function() {
  tibble::tribble(
    ~sample_id, ~subject, ~timepoint,
    "S01_T1", "S01", "T1", "S01_T2", "S01", "T2", "S01_T3", "S01", "T3",
    "S02_T1", "S02", "T1", "S02_T2", "S02", "T2", "S02_T3", "S02", "T3",
    "S28_T1", "S28", "T1",
    "S29_T2", "S29", "T2", "S29_T3", "S29", "T3"
  )
}

fixture_matrix <- function(meta) {
  matrix(
    seq_len(nrow(meta) * 2),
    nrow = 2, dimnames = list(c("f1", "f2"), meta$sample_id)
  )
}

test_that("subject_window subtracts in the stated direction", {
  source(here::here("functions", "association.R"))

  meta <- tibble::tribble(
    ~sample_id, ~subject, ~timepoint,
    "S01_T1", "S01", "T1", "S01_T2", "S01", "T2"
  )
  mat <- matrix(c(10, 100, 13, 106),
    nrow = 2, dimnames = list(c("f1", "f2"), meta$sample_id)
  )

  out <- subject_window(mat, "training", meta)
  expect_equal(colnames(out), "S01")
  expect_equal(unname(out[, "S01"]), c(3, 6))
})

test_that("subject_window drops subjects missing a needed timepoint", {
  source(here::here("functions", "association.R"))
  meta <- fixture_meta()
  mat <- fixture_matrix(meta)

  training <- subject_window(mat, "training", meta)
  expect_setequal(colnames(training), c("S01", "S02"))

  acute <- subject_window(mat, "acute", meta)
  expect_setequal(colnames(acute), c("S01", "S02", "S29"))
  expect_false("S28" %in% colnames(acute))

  total <- subject_window(mat, "total", meta)
  expect_setequal(colnames(total), c("S01", "S02"))
})

test_that("a level window returns the value itself, not a difference", {
  source(here::here("functions", "association.R"))
  meta <- fixture_meta()
  mat <- fixture_matrix(meta)

  # T1 exists for S01, S02 and S28 but not S29.
  t1 <- subject_window(mat, "T1", meta)
  expect_setequal(colnames(t1), c("S01", "S02", "S28"))
  expect_equal(
    unname(t1[, "S01"]),
    unname(mat[, "S01_T1"])
  )

  t3 <- subject_window(mat, "T3", meta)
  expect_setequal(colnames(t3), c("S01", "S02", "S29"))
  expect_false("S28" %in% colnames(t3))
})

test_that("total change equals training plus acute where both exist", {
  source(here::here("functions", "association.R"))
  meta <- fixture_meta()
  mat <- fixture_matrix(meta)

  both <- c("S01", "S02")
  tr <- subject_window(mat, "training", meta)[, both, drop = FALSE]
  ac <- subject_window(mat, "acute", meta)[, both, drop = FALSE]
  tot <- subject_window(mat, "total", meta)[, both, drop = FALSE]
  expect_equal(tot, tr + ac)
})

test_that("subject_window pairs each subject with its own timepoints", {
  source(here::here("functions", "association.R"))

  # Rows deliberately out of subject order, so a positional pairing would
  # silently subtract one subject's T1 from another's T2.
  meta <- tibble::tribble(
    ~sample_id, ~subject, ~timepoint,
    "S02_T2", "S02", "T2",
    "S01_T1", "S01", "T1",
    "S02_T1", "S02", "T1",
    "S01_T2", "S01", "T2"
  )
  mat <- matrix(c(20, 1, 30, 10),
    nrow = 1, dimnames = list("f1", meta$sample_id)
  )

  # S01 rises 1 -> 10, S02 falls 30 -> 20. A positional pairing would cross
  # them and return something else for both.
  out <- subject_window(mat, "training", meta)
  expect_equal(unname(out[, "S01"]), 9)
  expect_equal(unname(out[, "S02"]), -10)
})

test_that("subject_window rejects an unknown window", {
  source(here::here("functions", "association.R"))
  meta <- fixture_meta()
  expect_error(
    subject_window(fixture_matrix(meta), "nonsense", meta), "unknown window"
  )
})

test_that("associate recovers a slope it was given", {
  source(here::here("functions", "association.R"))

  set.seed(42)
  y <- setNames(rnorm(20), paste0("S", 1:20))
  feat <- rbind(
    flat = rnorm(20, sd = 0.1),
    sloped = 2 * y + rnorm(20, sd = 0.1)
  )
  colnames(feat) <- names(y)

  res <- associate(feat, y)
  expect_setequal(res$feature, c("flat", "sloped"))
  expect_equal(res$slope[res$feature == "sloped"], 2, tolerance = 0.1)
  expect_lt(res$bh[res$feature == "sloped"], 0.01)
  expect_gt(res$bh[res$feature == "flat"], 0.05)
  expect_equal(attr(res, "n"), 20)
})

test_that("associate drops subjects with a missing phenotype", {
  source(here::here("functions", "association.R"))

  y <- c(S1 = 1, S2 = 2, S3 = NA, S4 = 4, S5 = 5, S6 = 6, S7 = 7)
  feat <- matrix(rnorm(14),
    nrow = 2, dimnames = list(c("f1", "f2"), names(y))
  )

  res <- associate(feat, y)
  expect_equal(attr(res, "n"), 6)
})

test_that("associate aligns on names, not column order", {
  source(here::here("functions", "association.R"))

  set.seed(1)
  y <- setNames(seq_len(12), paste0("S", 1:12))
  feat <- rbind(sloped = 3 * unname(y))
  colnames(feat) <- names(y)

  shuffled <- feat[, rev(colnames(feat)), drop = FALSE]
  expect_equal(
    associate(feat, y)$slope,
    associate(shuffled, y)$slope,
    tolerance = 1e-10
  )
})
