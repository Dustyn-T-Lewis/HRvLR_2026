test_that("fit_hlm reproduces a direct lmerTest REML fit", {
  source(here::here("functions", "shared_hlm.R"))
  pacman::p_load(lmerTest, lme4)

  set.seed(1)
  d <- expand.grid(subject = factor(1:12), Timepoint = factor(c("T1", "T2")))
  d$Group <- factor(rep(c("HR", "LR"), each = 12))
  d$value <- rnorm(nrow(d)) + as.integer(d$Timepoint) * 0.5

  got <- fit_hlm(d, response = "value")
  ref <- lmerTest::lmer(
    value ~ Group * Timepoint + (1 | subject),
    data = d, REML = TRUE
  )
  expect_s4_class(got, "lmerModLmerTest")
  expect_equal(lme4::fixef(got), lme4::fixef(ref), tolerance = 1e-8)
  # Satterthwaite p-values survive (this is why lmerTest, not plain lme4)
  expect_true("Pr(>|t|)" %in% colnames(summary(got)$coefficients))
})
