test_that("storey_q keeps NA where p is NA and equals BH below 100 tests", {
  p <- c(0.001, 0.01, NA, 0.04, 0.2, 0.9)
  q <- storey_q(p)
  expect_length(q, length(p))
  expect_true(is.na(q[3]))
  expect_equal(q[-3], p.adjust(p[-3], "BH"))
})

test_that("storey_q matches qvalue's bootstrap estimate from 100 tests on", {
  set.seed(1)
  p <- c(runif(150), rbeta(50, 0.2, 5))
  expect_equal(storey_q(p), qvalue::qvalue(p, pi0.method = "bootstrap")$qvalues)
  expect_true(all(storey_q(p) <= p.adjust(p, "BH") + 1e-12))
})

test_that("storey_q falls back to BH when qvalue cannot fit", {
  local_mocked_bindings(qvalue = function(...) stop("cannot fit"), .package = "qvalue")
  set.seed(2)
  p <- runif(120)
  expect_equal(storey_q(p), p.adjust(p, "BH"))
})
