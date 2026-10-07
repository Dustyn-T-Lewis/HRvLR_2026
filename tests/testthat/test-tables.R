test_that("pairing_value looks a trait up by person", {
  pairings <- list(training = tibble(subject = c("A", "B"), pct_mcsa = c(5, 7)))
  expect_equal(unname(pairing_value(pairings, "training", "pct_mcsa", c("B", "A"))), c(7, 5))
})

test_that("status_of, cap and as_grid shape results for the overview maps", {
  expect_equal(
    status_of(c(0.001, 0.01, 0.5, NA), c(0.01, 0.2, 0.9, NA)),
    c("BH < 0.05", "p < 0.05", "not significant", NA)
  )
  expect_equal(cap(c(-100, 1:99)), unname(quantile(c(100, 1:99), 0.99)))
  d <- tibble(feature = c("a", "b"), contrast = c("x", "y"), effect = c(1, 2))
  grid <- as_grid(d, c("a", "b"), c("x", "y"), "contrast", \(d) d$effect)
  expect_equal(grid["a", "x"], 1)
  expect_true(is.na(grid["a", "y"]))
})

test_that("slots, capped_count and counts_line lay out report pages", {
  expect_length(slots(c("p1", "p2")), 12)
  expect_equal(capped_count(tibble(task = rep(c("a", "b"), c(30, 3))), "task"), 27)
  line <- counts_line(
    list(n_nominal = 3, n_tested = 10, n_fdr_05 = 1), "proteins", tibble(page = 1, n_pages = 2)
  )
  expect_equal(line, "3 of 10 proteins at p < 0.05, 1 at BH < 0.05. Page 1 of 2.")
})
