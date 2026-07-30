pacman::p_load(here, testthat, openxlsx, withr)
source(here("functions", "sweep_manifest.R"))

seed_cont_leaf <- function(root, level, config, phenotype, model, q2, p) {
  d <- file.path(root, level, config, phenotype, model, "c_data")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write_sweep_workbook(file.path(d, "results.xlsx"), list(
    summary = data.frame(
      level = level, config = config, outcome = phenotype, model = model,
      n = 15L, p = 79L, B = c(0L, 200L), q2 = q2,
      perm_p_q2 = c(NA, p), null_q2_mean = c(NA, -0.3), null_q2_sd = c(NA, 0.4)
    ),
    null = data.frame(outcome = phenotype, model = model, q2 = c(-0.2, -0.4)),
    predictions = data.frame(
      outcome = phenotype, model = model, subject = "HR_S03", y = 1, pred = 1
    ),
    selection = data.frame(
      outcome = phenotype, model = model,
      feature = c(
        "HALLMARK_MTORC1_SIGNALING@T1", "HALLMARK_MTORC1_SIGNALING@T2",
        "HALLMARK_COAGULATION"
      ),
      folds = 15L, freq = c(1, 1, 0.5)
    )
  ))
}

test_that("build_manifest reads the max-B row, never B = 0", {
  root <- withr::local_tempdir()
  seed_cont_leaf(root, "pathways", "T2", "d_mcsa", "spls", 0.41, 0.015)

  m <- build_manifest(root, "cont",
    screen_size = 792L, root_name = "F06_prediction"
  )

  expect_equal(nrow(m), 1L)
  expect_equal(m$metric_b200, 0.41, tolerance = 1e-8)
  expect_equal(m$perm_p, 0.015, tolerance = 1e-8)
  expect_equal(m$metric_b0, 0.41, tolerance = 1e-8)
  expect_equal(m$n_cells_screened, 792L)
})

test_that("build_manifest labels leakage by feature level", {
  root <- withr::local_tempdir()
  seed_cont_leaf(root, "pathways", "T2", "d_mcsa", "spls", 0.41, 0.015)
  seed_cont_leaf(root, "modules", "acute", "d_mcsa", "lasso", 0.55, 0.010)
  seed_cont_leaf(root, "proteins", "T1", "d_mcsa", "enet", -0.2, 0.7)

  m <- build_manifest(root, "cont",
    screen_size = 792L, root_name = "F06_prediction"
  )

  expect_equal(m$leakage[m$level == "pathways"], "leakage-free")
  expect_equal(m$leakage[m$level == "modules"], "optimistic")
  expect_equal(m$leakage[m$level == "proteins"], "optimistic")
})

test_that("build_manifest collapses trajectory suffixes in top_drivers", {
  root <- withr::local_tempdir()
  seed_cont_leaf(root, "pathways", "trajectory", "d_mcsa", "spls", 0.3, 0.04)

  m <- build_manifest(root, "cont",
    screen_size = 792L, root_name = "F06_prediction"
  )

  expect_equal(
    m$top_drivers,
    "HALLMARK_MTORC1_SIGNALING; HALLMARK_COAGULATION"
  )
  expect_equal(m$n_drivers, 2L)
})
