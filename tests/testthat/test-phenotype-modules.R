source(here::here(
  "03_Features", "05_phenotype_modules", "a_script", "trait_helpers.R"
))

test_that("cor_scan recovers a planted association", {
  set.seed(42)
  pheno <- rnorm(14)
  mat <- rbind(
    signal = pheno + rnorm(14, sd = 0.2),
    noise = rnorm(14)
  )
  scan <- cor_scan(mat, pheno, perm_index(14, 500, seed = 1))
  expect_gt(scan$obs$rho[1], 0.8)
  expect_lt(scan$obs$emp_p[1], 0.01)
  expect_gt(scan$obs$emp_p[2], 0.05)
  expect_equal(dim(scan$null_rho), c(2L, 500L))
})

test_that("cor_scan tolerates an NA in the phenotype", {
  set.seed(42)
  pheno <- c(rnorm(13), NA)
  mat <- matrix(rnorm(28), nrow = 2, dimnames = list(c("a", "b"), NULL))
  scan <- cor_scan(mat, pheno, perm_index(14, 100, seed = 1))
  expect_false(anyNA(scan$obs$rho))
  expect_false(anyNA(scan$obs$emp_p))
})

test_that("consistent_call needs sign agreement and two significant", {
  rho <- rbind(
    all_pos_two_sig = c(0.6, 0.7, 0.5),
    sign_flip = c(0.6, -0.7, 0.5),
    one_sig = c(0.6, 0.7, 0.5)
  )
  p <- rbind(
    all_pos_two_sig = c(0.01, 0.02, 0.30),
    sign_flip = c(0.01, 0.02, 0.03),
    one_sig = c(0.01, 0.40, 0.30)
  )
  expect_equal(
    unname(consistent_call(rho, p)),
    c(TRUE, FALSE, FALSE)
  )
})

test_that("null_emp_p gives the largest |rho| the smallest p", {
  null_rho <- matrix(c(0.9, 0.1, 0.5, -0.2), nrow = 1)
  p <- null_emp_p(null_rho)
  expect_equal(order(p[1, ]), order(-abs(null_rho[1, ])))
  expect_equal(min(p), 1 / 4, tolerance = 1e-8)
})

test_that("choose_k finds two planted clusters", {
  set.seed(42)
  base_a <- rnorm(20)
  base_b <- rnorm(20)
  x <- rbind(
    t(replicate(15, base_a + rnorm(20, sd = 0.3))),
    t(replicate(15, base_b + rnorm(20, sd = 0.3)))
  )
  d <- as.dist(1 - cor(t(x)))
  expect_equal(choose_k(d, hclust(d, "ward.D2"), 2:6), 2L)
})
