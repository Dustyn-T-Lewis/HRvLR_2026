pacman::p_load(here, testthat, openxlsx, withr)
source(here("functions", "sweep_split.R"))

fake_pred_workbook <- function(dir) {
  path <- file.path(dir, "results.xlsx")
  write_sweep_workbook(path, list(
    summary = data.frame(
      level = "pathways", config = "T2",
      outcome = c("d_mcsa", "d_mcsa", "d_fcsa_I"),
      model = "spls", n = 15L, p = 79L, B = c(0L, 200L, 200L),
      q2 = c(0.41, 0.41, -0.09), perm_p_q2 = c(NA, 0.015, 0.194)
    ),
    null = data.frame(
      outcome = c("d_mcsa", "d_fcsa_I"), model = "spls", q2 = c(-0.2, -0.3)
    ),
    predictions = data.frame(
      outcome = c("d_mcsa", "d_fcsa_I"), model = "spls",
      subject = "HR_S03", y = c(1.2, 3.4), pred = c(1.0, 3.0)
    ),
    selection = data.frame(
      outcome = c("d_mcsa", "d_fcsa_I"), model = "spls",
      feature = c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_COAGULATION"),
      folds = 15L, freq = 1
    )
  ))
  path
}

test_that("split_pred_leaf writes one leaf per outcome and preserves values", {
  root_dir <- withr::local_tempdir()
  src <- withr::local_tempdir()
  path <- fake_pred_workbook(src)

  written <- split_pred_leaf(path, root_dir, "pathways", "T2", "spls")

  expect_setequal(written, c("d_mcsa", "d_fcsa_I"))

  out <- file.path(
    root_dir, "pathways", "T2", "d_mcsa", "spls",
    "c_data", "results.xlsx"
  )
  expect_true(file.exists(out))
  expect_setequal(
    openxlsx::getSheetNames(out),
    c("summary", "null", "predictions", "selection")
  )

  s <- openxlsx::read.xlsx(out, "summary")
  expect_equal(nrow(s), 2L)
  expect_true(all(s$outcome == "d_mcsa"))
  expect_equal(s$q2[s$B == 200], 0.41, tolerance = 1e-8)
  expect_equal(s$perm_p_q2[s$B == 200], 0.015, tolerance = 1e-8)

  sel <- openxlsx::read.xlsx(out, "selection")
  expect_equal(sel$feature, "HALLMARK_MTORC1_SIGNALING")
})

test_that(paste(
  "split_pred_leaf with an explicit phenotype relocates a",
  "single-outcome leaf"
), {
  root_dir <- withr::local_tempdir()
  src <- withr::local_tempdir()
  path <- file.path(src, "results.xlsx")
  write_sweep_workbook(path, list(
    summary = data.frame(
      level = "modules", config = "total", model = "plain", metric = "AUC",
      n = 15L, p = 13L, B = c(0L, 200L), estimate = 0.78,
      perm_p = c(NA, 0.045)
    ),
    null = data.frame(model = "plain", auc = c(0.4, 0.6)),
    predictions = data.frame(
      model = "plain", subject = "HR_S03", y = 1L, pred = 0.7
    ),
    selection = data.frame(
      model = "plain", feature = "ME_green", folds = 15L, freq = 1
    )
  ))

  written <- split_pred_leaf(path, root_dir, "modules", "total", "plain",
    phenotype = "HR_LR"
  )

  expect_equal(written, "HR_LR")
  out <- file.path(
    root_dir, "modules", "total", "HR_LR", "plain",
    "c_data", "results.xlsx"
  )
  s <- openxlsx::read.xlsx(out, "summary")
  expect_equal(s$estimate[s$B == 200], 0.78, tolerance = 1e-8)
})

test_that("split_assoc_leaf turns each outcome sheet into a leaf", {
  root_dir <- withr::local_tempdir()
  src <- withr::local_tempdir()
  path <- file.path(src, "results.xlsx")
  cell <- function(oc) {
    data.frame(
      outcome = oc, n = 15L,
      feature = c("HALLMARK_ADIPOGENESIS", "HALLMARK_COAGULATION"),
      effect = c(0.001, -0.002), t = c(1.02, -1.81), p = c(0.32, 0.09),
      bh = c(0.46, 0.22)
    )
  }
  write_sweep_workbook(path, list(
    d_mcsa = cell("d_mcsa"), group_diff = cell("group_diff"),
    summary = data.frame(
      level = "pathways", config = "T2", method = "limma",
      outcome = c("d_mcsa", "group_diff"), n_feature = 2L, n_nominal = 0L,
      n_bh = 0L, med_abs_effect = 0.0015
    )
  ))

  written <- split_assoc_leaf(path, root_dir, "pathways", "T2", "limma")

  expect_setequal(written, c("d_mcsa", "group_diff"))
  out <- file.path(
    root_dir, "pathways", "T2", "group_diff", "limma",
    "c_data", "results.xlsx"
  )
  expect_setequal(openxlsx::getSheetNames(out), c("cell", "summary"))
  d <- openxlsx::read.xlsx(out, "cell")
  expect_equal(nrow(d), 2L)
  expect_equal(d$t, c(1.02, -1.81), tolerance = 1e-8)
  expect_equal(openxlsx::read.xlsx(out, "summary")$outcome, "group_diff")
})

test_that("split_pred_leaf errors when the workbook has no usable outcome", {
  root_dir <- withr::local_tempdir()
  src <- withr::local_tempdir()
  path <- file.path(src, "results.xlsx")
  write_sweep_workbook(path, list(
    summary = data.frame(
      level = "pathways", config = "T2", model = "spls", B = 200L
    )
  ))
  expect_error(
    split_pred_leaf(path, root_dir, "pathways", "T2", "spls"),
    "no outcome"
  )
})
