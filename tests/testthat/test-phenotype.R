test_that("phenotype_tables computes % change and pairs each window with its trait table", {
  traits <- c("fcsa_mixed", "fcsa_type1", "fcsa_type2", "mcsa", "leg_press_1rm", "leg_ext_1rm")
  phenotype <- tibble(subject = c("HR_S1", "LR_S2"), arm = c("HR", "LR"))
  for (t in traits) {
    phenotype[[paste0(t, "_pre_u")]] <- c(100, 200)
    phenotype[[paste0(t, "_post_u")]] <- c(110, 150)
  }
  phenotype$comp_hypertrophy <- c(1, -1)
  phenotype$volume_load_total_kg <- c(5000, 6000)
  out <- phenotype_tables(phenotype)
  expect_equal(out$pairings$training$pct_fcsa_mixed, c(10, -25))
  expect_identical(out$pairings$training, out$pairings$acute)
  expect_equal(out$pairings$baseline_level$mcsa, c(100, 200))
  expect_equal(out$pairings$post_bout_level$mcsa, c(110, 150))
  expect_named(out$pairings, c(
    "training", "baseline", "acute", "baseline_level", "post_level", "post_bout_level"
  ))
})
