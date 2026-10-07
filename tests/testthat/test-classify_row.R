spec_of <- function(pos, neg) {
  list(
    positive = matrix(pos, 1, dimnames = list("f")),
    negative = matrix(neg, 1, dimnames = list("f"))
  )
}

test_that("classify_row gives the Mann-Whitney AUC, Wilcoxon p and DeLong interval", {
  pos <- c(5.1, 6.3, 4.8, 7.0, 5.9)
  neg <- c(4.2, 5.0, 3.9, 6.1, 4.4, 4.0)
  out <- classify_row(spec_of(pos, neg), "f")
  test <- wilcox.test(pos, neg)
  delong <- pROC::roc(controls = neg, cases = pos, direction = "<", quiet = TRUE) |>
    pROC::ci.auc(method = "delong") |>
    as.numeric()
  expect_equal(unname(out["auc"]), unname(test$statistic) / (5 * 6))
  expect_equal(unname(out["p"]), test$p.value)
  expect_equal(unname(out[c("auc_lower", "auc_upper")]), delong[c(1, 3)])
  expect_equal(unname(out[c("n_positive", "n_negative")]), c(5, 6))
})

test_that("classify_row keeps an AUC below 0.5 instead of flipping it", {
  out <- classify_row(spec_of(c(1, 2, 3, 4.5), c(4, 5, 6, 7)), "f")
  expect_lt(out[["auc"]], 0.5)
  expect_lte(out[["auc_lower"]], out[["auc"]])
})

test_that("classify_row gives no interval at perfect separation and nothing below three a side", {
  perfect <- classify_row(spec_of(c(10, 11, 12), c(1, 2, 3)), "f")
  expect_equal(perfect[["auc"]], 1)
  expect_true(is.na(perfect[["auc_lower"]]) && is.na(perfect[["auc_upper"]]))
  thin <- classify_row(spec_of(c(1, 2, NA), c(3, 4, 5)), "f")
  expect_true(is.na(thin[["auc"]]) && is.na(thin[["p"]]))
  expect_equal(thin[["n_positive"]], 2)
})
