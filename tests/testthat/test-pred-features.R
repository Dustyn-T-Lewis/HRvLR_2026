test_that("pred_contrast_matrix drops subjects missing a required timepoint", {
  source(here::here("functions", "pred_features.R"))

  # S01, S02 have all three timepoints; S28-style has T1 only; S29-style has
  # T2/T3 only -- the real roster's own pattern.
  meta <- tibble::tribble(
    ~sample, ~subject, ~timepoint,
    "S01_T1", "S01", "T1", "S01_T2", "S01", "T2", "S01_T3", "S01", "T3",
    "S02_T1", "S02", "T1", "S02_T2", "S02", "T2", "S02_T3", "S02", "T3",
    "S28_T1", "S28", "T1",
    "S29_T2", "S29", "T2", "S29_T3", "S29", "T3"
  )
  feature_mat <- matrix(
    seq_len(nrow(meta) * 2),
    nrow = 2, dimnames = list(c("f1", "f2"), meta$sample)
  )

  baseline <- pred_contrast_matrix(feature_mat, meta, "T1")
  expect_setequal(rownames(baseline), c("S01", "S02", "S28"))

  training <- pred_contrast_matrix(feature_mat, meta, "training")
  expect_setequal(rownames(training), c("S01", "S02"))
  expect_false("S28" %in% rownames(training))
  expect_false("S29" %in% rownames(training))

  acute <- pred_contrast_matrix(feature_mat, meta, "acute")
  expect_setequal(rownames(acute), c("S01", "S02", "S29"))
  expect_false("S28" %in% rownames(acute))

  trajectory <- pred_contrast_matrix(feature_mat, meta, "trajectory")
  expect_setequal(rownames(trajectory), c("S01", "S02"))
  expect_equal(ncol(trajectory), 2 * 3)
})

test_that("pred_contrast_matrix computes the diff in the stated direction", {
  source(here::here("functions", "pred_features.R"))

  meta <- tibble::tribble(
    ~sample, ~subject, ~timepoint,
    "S01_T1", "S01", "T1", "S01_T2", "S01", "T2"
  )
  feature_mat <- matrix(c(10, 100, 13, 106),
    nrow = 2, dimnames = list(c("f1", "f2"), meta$sample)
  )

  out <- pred_contrast_matrix(feature_mat, meta, "training")
  expect_equal(unname(out["S01", ]), c(3, 6))
})

test_that("pred_contrast_matrix rejects an unknown contrast", {
  source(here::here("functions", "pred_features.R"))
  meta <- tibble::tibble(sample = "S01_T1", subject = "S01", timepoint = "T1")
  feature_mat <- matrix(1, dimnames = list("f1", "S01_T1"))
  expect_error(pred_contrast_matrix(feature_mat, meta, "nonsense"))
})

test_that("pred_eigengene_matrix reshapes long eigengenes to wide correctly", {
  source(here::here("functions", "pred_features.R"))

  long <- tibble::tribble(
    ~sample_id, ~group_id, ~ME,
    "S01_T1", "blue", 0.1, "S01_T1", "red", -0.2,
    "S02_T1", "blue", 0.4, "S02_T1", "red", 0.5
  )
  tmp <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(long, tmp)

  wide <- pred_eigengene_matrix(tmp, c("S01_T1", "S02_T1"))
  expect_equal(dim(wide), c(2, 2))
  expect_equal(rownames(wide), c("ME_blue", "ME_red"))
  expect_equal(wide["ME_blue", "S01_T1"], 0.1)
  expect_equal(wide["ME_red", "S02_T1"], 0.5)
})

test_that("pred_outcome encodes group as 0/1 and reads phenotype columns", {
  source(here::here("functions", "pred_features.R"))

  bundle <- list(
    meta = tibble::tribble(
      ~subject, ~group,
      "S01", "HR", "S02", "LR"
    ),
    pheno = tibble::tribble(
      ~subject, ~d_mcsa,
      "S01", 120, "S02", -30
    )
  )

  grp <- pred_outcome(bundle, "group")
  expect_equal(unname(grp[["S01"]]), 1L)
  expect_equal(unname(grp[["S02"]]), 0L)

  mcsa <- pred_outcome(bundle, "d_mcsa")
  expect_equal(unname(mcsa[["S01"]]), 120)
})

test_that("pred_paths points each tree at its own eigengene file", {
  source(here::here("functions", "pred_features.R"))

  cat_path <- pred_paths("categorical")$eigen
  cont_path <- pred_paths("continuous")$eigen
  expect_match(cat_path, "categorical")
  expect_match(cont_path, "continuous")
  expect_false(cat_path == cont_path)
})
