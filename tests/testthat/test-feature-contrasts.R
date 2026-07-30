test_that("fit_feature_contrasts reproduces the stage 03 proteoDA fit", {
  source(here::here("functions", "feature_contrasts.R"))

  # The 34 proteins the missingness filter admits but the model cannot test are
  # documented in 03_DEP/contrasts.R; limma warns about them by design.
  expect_warning(
    got <- verify_protein_equivalence(),
    "Partial NA coefficients for 34 probe"
  )

  expect_equal(nrow(got), 9L)
  expect_true(all(got$equivalent))
  expect_true(all(got$n == 1900L))
  expect_lt(max(got$max_d_logFC), 1e-10)
  expect_lt(max(got$max_d_p), 1e-10)
  expect_lt(max(got$max_d_bh), 1e-10)
})

test_that("fit_feature_contrasts returns nine contrasts, BH within each", {
  source(here::here("functions", "feature_contrasts.R"))
  pacman::p_load(dplyr)

  meta <- feature_metadata()
  set.seed(1)
  mat <- matrix(rnorm(20 * nrow(meta)),
    nrow = 20, dimnames = list(paste0("f", 1:20), meta$sample_id)
  )

  got <- fit_feature_contrasts(mat, meta)

  expect_setequal(
    unique(got$contrast), trimws(sub("=.*$", "", HRVLR_CONTRASTS))
  )
  expect_equal(nrow(got), 9L * 20L)
  expect_true(all(got$bh >= got$p - 1e-12))

  # BH within a contrast never sees the other eight: each contrast's smallest q
  # is its own smallest p times its own feature count, not 180.
  per <- got |>
    group_by(.data$contrast) |>
    summarise(
      ratio = min(.data$bh) / min(.data$p), .groups = "drop"
    )
  expect_true(all(per$ratio <= 20 + 1e-8))
})

test_that("fit_feature_contrasts rejects samples the design cannot place", {
  source(here::here("functions", "feature_contrasts.R"))

  meta <- feature_metadata()
  mat <- matrix(0,
    nrow = 2, ncol = 2,
    dimnames = list(c("f1", "f2"), c(meta$sample_id[1], "S99_T9"))
  )

  expect_error(
    fit_feature_contrasts(mat, meta), "absent from the normalisation"
  )
})
