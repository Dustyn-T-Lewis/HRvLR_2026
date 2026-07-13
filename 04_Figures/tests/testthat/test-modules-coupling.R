pacman::p_load(here, testthat, dplyr)
source(here("04_Figures", "modules", "a_script", "functions", "coupling.R"))

test_that("coupling recovers a planted correlation and BH-adjusts", {
  subj <- paste0("S", 1:12)
  me1 <- rnorm(12)
  pheno_val <- me1 * 2 + rnorm(12, sd = 0.1)
  me <- data.frame(
    MEred = me1, MEblue = rnorm(12),
    row.names = paste0(subj, "_T1")
  )
  meta <- data.frame(sample_id = paste0(subj, "_T1"), subject = subj, timepoint = "T1")
  pheno <- data.frame(subject = subj, d_mcsa = pheno_val, comp_hypertrophy = rnorm(12))

  out <- module_coupling(me, meta, pheno, "baseline")
  red <- dplyr::filter(out, module == "MEred", trait == "d_mcsa")
  expect_gt(red$r, 0.8)
  expect_true(all(c("p_bh", "phase") %in% names(out)))
  expect_equal(unique(out$phase), "baseline")
})
