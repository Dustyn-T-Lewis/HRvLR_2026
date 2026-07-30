pacman::p_load(here, testthat, dplyr, limma)
source(here("functions", "feature_contrasts.R"))

synthetic <- function(seed = 31) {
  set.seed(seed)
  subjects <- sprintf("S%02d", 1:16)
  meta <- expand.grid(
    subject = subjects, tp = c("T1", "T2", "T3"), stringsAsFactors = FALSE
  ) |>
    mutate(
      arm = ifelse(match(.data$subject, subjects) <= 8, "HR", "LR"),
      group = factor(paste(.data$arm, .data$tp, sep = "_"),
        levels = GROUP_LEVELS
      ),
      sample_id = paste(.data$subject, .data$tp, sep = "_")
    ) |>
    select("sample_id", "group", "subject")
  mat <- matrix(rnorm(30 * nrow(meta)), nrow = 30)
  rownames(mat) <- paste0("p", 1:30)
  colnames(mat) <- meta$sample_id
  list(meta = meta, mat = mat)
}

test_that("an adjusted design keeps the six group columns and adds the covariate", {
  d <- synthetic()
  covar <- setNames(rnorm(ncol(d$mat)), colnames(d$mat))
  parts <- feature_design(d$mat, d$meta, adjust = covar)

  expect_true(all(GROUP_LEVELS %in% colnames(parts$design)))
  expect_true("adjust" %in% colnames(parts$design))
  expect_identical(ncol(parts$design), 7L)
})

test_that("the nine contrasts put no weight on the covariate", {
  d <- synthetic()
  covar <- setNames(rnorm(ncol(d$mat)), colnames(d$mat))
  parts <- feature_design(d$mat, d$meta, adjust = covar)

  expect_identical(
    colnames(parts$contrasts), trimws(sub("=.*$", "", HRVLR_CONTRASTS))
  )
  expect_true(all(parts$contrasts["adjust", ] == 0))
})

test_that("adjusting removes a signal the covariate fully explains", {
  d <- synthetic()
  covar <- setNames(rnorm(ncol(d$mat), mean = 5), colnames(d$mat))
  # Ten proteins are pure covariate plus noise, with no group structure at all.
  d$mat[1:10, ] <- rep(covar, each = 10) * 2 +
    matrix(rnorm(10 * ncol(d$mat), sd = 0.2), nrow = 10)

  plain <- fit_feature_contrasts(d$mat, d$meta)
  adj <- fit_feature_contrasts(d$mat, d$meta, adjust = covar)

  driven <- paste0("p", 1:10)
  expect_lt(
    median(adj$p[adj$feature %in% driven]),
    1
  )
  expect_gt(
    mean(abs(adj$logFC[adj$feature %in% driven])),
    0
  )
  expect_lt(
    mean(abs(adj$logFC[adj$feature %in% driven])),
    mean(abs(plain$logFC[plain$feature %in% driven]))
  )
})

test_that("adjusting by a covariate orthogonal to the design changes little", {
  d <- synthetic()
  covar <- setNames(rnorm(ncol(d$mat)), colnames(d$mat))
  plain <- fit_feature_contrasts(d$mat, d$meta)
  adj <- fit_feature_contrasts(d$mat, d$meta, adjust = covar)

  expect_equal(plain$logFC, adj$logFC, tolerance = 0.35)
  expect_setequal(plain$feature, adj$feature)
  expect_identical(nrow(plain), nrow(adj))
})

test_that("a covariate missing a sample is rejected rather than recycled", {
  d <- synthetic()
  covar <- setNames(rnorm(ncol(d$mat)), colnames(d$mat))
  expect_error(
    feature_design(d$mat, d$meta, adjust = covar[-1]),
    "adjust"
  )
})
