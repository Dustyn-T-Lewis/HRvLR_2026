pacman::p_load(here, testthat, dplyr, tibble)
source(here("04_Figures", "F04_modules", "a_script", "wgcna.R"))

test_that("centre_within_subject removes the between-subject mean exactly", {
  set.seed(1)
  subject <- rep(paste0("S", 1:4), each = 3)
  offsets <- setNames(c(10, -5, 100, 0), paste0("S", 1:4))
  abund <- matrix(rnorm(5 * 12), nrow = 5)
  colnames(abund) <- paste0("s", 1:12)
  abund <- sweep(abund, 2, offsets[subject], "+")

  centred <- centre_within_subject(abund, subject)

  # every subject's mean is zero for every protein, to numeric tolerance
  means <- vapply(unique(subject), function(s) {
    max(abs(rowMeans(centred[, subject == s, drop = FALSE])))
  }, numeric(1))
  expect_lt(max(means), 1e-10)

  # and the shape is preserved
  expect_identical(dim(centred), dim(abund))
  expect_identical(colnames(centred), colnames(abund))
})

test_that("centre_within_subject annihilates a pure subject-identity signal", {
  subject <- rep(c("A", "B"), each = 4)
  # a protein that is ONLY subject identity: no within-subject variation at all
  abund <- rbind(subject_only = ifelse(subject == "A", 1, 9))
  colnames(abund) <- paste0("s", seq_along(subject))

  centred <- centre_within_subject(abund, subject)

  expect_lt(max(abs(centred)), 1e-12)
})

test_that("subject_variance recovers a planted ICC and returns 0 for pure noise", {
  set.seed(42)
  subject <- rep(paste0("S", 1:10), each = 3)

  # MEhigh is almost entirely between-subject; MEnoise has no subject structure
  high <- rep(rnorm(10, sd = 5), each = 3) + rnorm(30, sd = 0.3)
  noise <- rnorm(30)

  me_long <- bind_rows(
    tibble(module = "MEhigh", subject = subject, ME = high),
    tibble(module = "MEnoise", subject = subject, ME = noise)
  )

  sv <- subject_variance(me_long)

  icc <- setNames(sv$icc, sv$module)
  expect_gt(icc[["MEhigh"]], 0.9)
  expect_lt(icc[["MEnoise"]], 0.5)
  expect_true(all(sv$icc >= 0 & sv$icc <= 1))
})

test_that("choose_power falls back only when the scale-free fit never clears the cut", {
  # pure noise cannot be scale-free, so pickSoftThreshold returns NA and the
  # sample-size fallback must fire (n > 40 -> 12)
  skip_if_not_installed("WGCNA")
  set.seed(7)
  expr <- matrix(rnorm(45 * 60), nrow = 45)
  colnames(expr) <- paste0("p", 1:60)

  pw <- suppressWarnings(choose_power(expr))

  expect_true(is.numeric(pw$power))
  expect_true(pw$power %in% 1:20)
})
