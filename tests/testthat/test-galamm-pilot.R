source(here::here(
  "03_Features", "04_galamm_pilot", "a_script", "pilot_helpers.R"
))

test_that("varident_ratios recovers simulated per-timepoint SDs", {
  set.seed(42)
  n_subj <- 1000
  d <- expand.grid(
    subject = factor(seq_len(n_subj)),
    timepoint = factor(c("T1", "T2", "T3"))
  )
  sds <- c(T1 = 1, T2 = 2, T3 = 0.5)
  d$y <- rnorm(n_subj)[d$subject] +
    rnorm(nrow(d), sd = sds[as.character(d$timepoint)])
  fit <- nlme::lme(y ~ timepoint,
    random = ~ 1 | subject, data = d,
    weights = nlme::varIdent(form = ~ 1 | timepoint)
  )
  ratios <- varident_ratios(fit)
  expect_named(ratios, c("T2", "T3"))
  expect_equal(unname(ratios["T2"]), 2, tolerance = 0.15)
  expect_equal(unname(ratios["T3"]), 0.5, tolerance = 0.15)

  rel <- coef(fit$modelStruct$varStruct,
    unconstrained = FALSE, allCoef = TRUE
  )
  expect_equal(unname(ratios["T2"]), unname(rel["T2"] / rel["T1"]),
    tolerance = 1e-8
  )
})

test_that("variance_explained matches the factor-model identity", {
  expect_equal(variance_explained(1, 1, 1), 0.5, tolerance = 1e-8)
  expect_equal(variance_explained(0, 1, 1), 0, tolerance = 1e-8)
  expect_equal(variance_explained(2, 0.5, 1), 2 / 3, tolerance = 1e-8)
})

test_that("emp_p applies the add-one permutation correction", {
  expect_equal(emp_p(0, 200), 1 / 201, tolerance = 1e-8)
  expect_equal(emp_p(200, 200), 1, tolerance = 1e-8)
})
