test_that("spearman_by_row matches cor.test and DescTools' Fisher-z interval", {
  set.seed(3)
  x <- matrix(rnorm(24), 2, 12, dimnames = list(c("a", "b"), NULL))
  outcome <- rnorm(12)
  outcome[4] <- NA
  out <- spearman_by_row(x, outcome)
  ref <- cor.test(x[1, -4], outcome[-4], method = "spearman", exact = FALSE)
  expect_equal(out$feature, c("a", "b"))
  expect_equal(out$n_subjects, c(11, 11))
  expect_equal(out$rho[1], unname(ref$estimate))
  expect_equal(out$p[1], ref$p.value)
  expect_equal(
    c(out$rho_lower[1], out$rho_upper[1]), unname(DescTools::CorCI(out$rho[1], 11)[2:3])
  )
  expect_true(all(out$rho_lower <= out$rho & out$rho <= out$rho_upper))
})

test_that("spearman_by_row gives no interval below four people", {
  out <- spearman_by_row(matrix(c(1, 2, 3), 1, dimnames = list("a", NULL)), c(3, 1, 2))
  expect_true(is.na(out$rho_lower) && is.na(out$rho_upper))
})

test_that("subject_fit returns the least-squares slope with an interval around it", {
  set.seed(4)
  predictors <- data.frame(trait = c(rnorm(9), NA))
  x <- matrix(rnorm(30), 3, 10, dimnames = list(c("a", "b", "c"), NULL))
  x[1, ] <- 2 * predictors$trait + rnorm(10, sd = 0.1)
  out <- subject_fit(x, predictors, ~trait, "trait")
  ols <- coef(lm(x[1, 1:9] ~ predictors$trait[1:9]))[[2]]
  expect_equal(out$estimate[out$feature == "a"], ols)
  expect_equal(unique(out$n_subjects), 9)
  expect_true(all(out$estimate_lower <= out$estimate & out$estimate <= out$estimate_upper))
  expect_lt(out$p[out$feature == "a"], 0.001)
})
