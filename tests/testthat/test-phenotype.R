pacman::p_load(here, testthat, readr)
source(here("00_input", "a_script", "build_phenotype.R"))

test_that("phenotype table has five traits and drops MyoVision", {
  tbl <- build_phenotype_table(here("00_input", "HRvLR_meta.csv"))
  expect_false("d_myovision_fcsa_I" %in% names(tbl))
  expect_true(all(c(
    "subject", "group_arm", "comp_hypertrophy",
    "d_fcsa_I", "d_fcsa_II", "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
  ) %in% names(tbl)))
  expect_equal(nrow(tbl), length(unique(tbl$subject)))
})
