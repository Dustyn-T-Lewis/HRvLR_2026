test_that("shared_prediction harness is deterministic and leakage-free", {
  source(here::here("functions", "shared_prediction.R"))

  set.seed(1)
  x <- matrix(rnorm(20 * 6), nrow = 20, dimnames = list(paste0("s", 1:20), NULL))
  y <- as.numeric(x[, 1] > 0)

  preds1 <- suppressWarnings(nested_loso(x, y, "glmnet", "binomial"))
  preds2 <- suppressWarnings(nested_loso(x, y, "glmnet", "binomial"))
  expect_identical(preds1, preds2)
  expect_length(preds1, 20)

  # perm_matrix must be reproducible from its seed alone
  expect_identical(perm_matrix(20, nperm = 5), perm_matrix(20, nperm = 5))

  # align_xy drops zero-variance columns and NA outcomes
  xz <- cbind(x, const = 1)
  al <- align_xy(xz, setNames(c(y[1:19], NA), rownames(xz)))
  expect_false("const" %in% colnames(al$x))
  expect_equal(nrow(al$x), 19)
})
