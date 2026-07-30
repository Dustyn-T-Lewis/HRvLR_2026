pacman::p_load(here, testthat, dplyr)
source(here("functions", "pi_permutation.R"))

# 16 subjects, 8 per arm, three timepoints each - the roster shape without the
# two partial subjects, which the permutation handles the same way.
synthetic_meta <- function() {
  subjects <- sprintf("S%02d", 1:16)
  arms <- rep(c("HR", "LR"), each = 8)
  expand.grid(
    subject = subjects, timepoint = c("T1", "T2", "T3"),
    stringsAsFactors = FALSE
  ) |>
    mutate(
      arm = arms[match(.data$subject, subjects)],
      group = factor(paste(.data$arm, .data$timepoint, sep = "_"),
        levels = GROUP_LEVELS
      ),
      sample_id = paste(.data$subject, .data$timepoint, sep = "_")
    ) |>
    select("sample_id", "group", "subject") |>
    arrange(.data$sample_id)
}

synthetic_matrix <- function(meta, effect = 0, n = 40, seed = 11) {
  set.seed(seed)
  mat <- matrix(rnorm(n * nrow(meta)), nrow = n)
  rownames(mat) <- paste0("p", seq_len(n))
  colnames(mat) <- meta$sample_id
  if (effect != 0) {
    is_hr <- grepl("^HR", meta$group)
    mat[1:10, is_hr] <- mat[1:10, is_hr] + effect
  }
  mat
}

test_that("pi_score is Xiao Eq.2 and pivots at two-fold", {
  expect_equal(pi_score(0.05, 1), 0.05)
  expect_equal(pi_score(0.05, 2), 0.0025)
  expect_gt(pi_score(0.05, 0.5), 0.05)
  expect_lt(pi_score(0.05, 3), 0.05)
  expect_equal(pi_score(0.05, -2), 0.0025)
})

test_that("pi_score propagates NA rather than inventing a value", {
  expect_true(is.na(pi_score(NA_real_, 2)))
  expect_true(is.na(pi_score(0.05, NA_real_)))
  expect_equal(pi_score(c(0.05, NA), c(1, 1)), c(0.05, NA))
})

test_that("permute_arm_labels keeps one arm per subject and the 8/8 balance", {
  meta <- synthetic_meta()
  set.seed(3)
  perm <- permute_arm_labels(meta)

  arm <- sub("_.*$", "", perm$group)
  per_subject <- tapply(arm, perm$subject, function(x) length(unique(x)))
  expect_true(all(per_subject == 1L))
  expect_equal(sort(table(unique(data.frame(perm$subject, arm))$arm)),
    sort(table(c(rep("HR", 8), rep("LR", 8)))),
    ignore_attr = TRUE
  )
})

test_that("permute_arm_labels leaves sample order and timepoints untouched", {
  meta <- synthetic_meta()
  set.seed(5)
  perm <- permute_arm_labels(meta)

  expect_identical(perm$sample_id, meta$sample_id)
  expect_identical(perm$subject, meta$subject)
  expect_identical(
    sub("^[^_]*_", "", as.character(perm$group)),
    sub("^[^_]*_", "", as.character(meta$group))
  )
  expect_identical(levels(perm$group), GROUP_LEVELS)
})

test_that("a planted arm effect beats its own permutation null", {
  meta <- synthetic_meta()
  mat <- synthetic_matrix(meta, effect = 3)
  got <- pi_permutation_null(mat, meta, n_perm = 30, seed = 1)

  base <- got[got$contrast == "Baseline_HRvLR", ]
  expect_gt(base$observed, base$null_median)
  expect_lt(base$emp_p, 0.05)
})

test_that("pure noise does not beat its own permutation null", {
  meta <- synthetic_meta()
  mat <- synthetic_matrix(meta, effect = 0)
  got <- pi_permutation_null(mat, meta, n_perm = 30, seed = 1)

  base <- got[got$contrast == "Baseline_HRvLR", ]
  expect_gt(base$emp_p, 0.05)
})

test_that("pi_permutation_null reports every contrast with a bounded emp_p", {
  meta <- synthetic_meta()
  mat <- synthetic_matrix(meta)
  got <- pi_permutation_null(mat, meta, n_perm = 10, seed = 2)

  expect_setequal(got$contrast, trimws(sub("=.*$", "", HRVLR_CONTRASTS)))
  expect_true(all(got$emp_p > 0 & got$emp_p <= 1))
  expect_true(all(got$n_perm == 10L))
})

test_that("pi_permutation_null is reproducible from its seed", {
  meta <- synthetic_meta()
  mat <- synthetic_matrix(meta)
  a <- pi_permutation_null(mat, meta, n_perm = 10, seed = 99)
  b <- pi_permutation_null(mat, meta, n_perm = 10, seed = 99)
  expect_equal(a, b)
})

test_that("the within-arm contrasts are permuted too, not held fixed", {
  meta <- synthetic_meta()
  mat <- synthetic_matrix(meta)
  got <- pi_permutation_null(mat, meta, n_perm = 10, seed = 4)
  expect_true(all(c("Training_HR", "Acute_LR") %in% got$contrast))
  expect_true(all(got$null_median[got$contrast == "Training_HR"] >= 0))
})
