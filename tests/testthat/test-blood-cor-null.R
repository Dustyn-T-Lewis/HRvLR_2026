pacman::p_load(here, testthat)
source(here("functions", "blood_cor_null.R"))

# A matrix where nothing tracks the index, plus two planted trackers. The null
# quantile has to sit above the noise and below the plants.
planted <- function(n_prot = 200, n_samp = 45, seed = 13) {
  set.seed(seed)
  mat <- matrix(rnorm(n_prot * n_samp), nrow = n_prot)
  rownames(mat) <- paste0("p", seq_len(n_prot))
  index <- rnorm(n_samp, mean = 27, sd = 1)
  mat[1, ] <- index * 2 + rnorm(n_samp, sd = 0.1)
  mat[2, ] <- -index * 2 + rnorm(n_samp, sd = 0.1)
  list(mat = mat, index = index)
}

test_that("blood_cor_null_quantile returns one finite correlation", {
  d <- planted()
  got <- blood_cor_null_quantile(d$mat, d$index, n_perm = 30)

  expect_length(got, 1L)
  expect_true(is.finite(got))
  expect_gt(got, 0)
  expect_lt(got, 1)
})

test_that("the null sits above unrelated proteins and below planted trackers", {
  d <- planted()
  cut <- blood_cor_null_quantile(d$mat, d$index, n_perm = 60)
  observed <- abs(as.vector(
    stats::cor(t(d$mat), d$index, method = "spearman")
  ))

  expect_gt(observed[1], cut)
  expect_gt(observed[2], cut)
  expect_lt(mean(observed[-(1:2)]), cut)
})

test_that("a higher probability gives a higher cut", {
  d <- planted()
  expect_gt(
    blood_cor_null_quantile(d$mat, d$index, n_perm = 40, prob = 0.999),
    blood_cor_null_quantile(d$mat, d$index, n_perm = 40, prob = 0.95)
  )
})

test_that("blood_cor_null_quantile is reproducible from its seed", {
  d <- planted()
  a <- blood_cor_null_quantile(d$mat, d$index, n_perm = 25, seed = 7)
  b <- blood_cor_null_quantile(d$mat, d$index, n_perm = 25, seed = 7)
  c <- blood_cor_null_quantile(d$mat, d$index, n_perm = 25, seed = 8)

  expect_identical(a, b)
  expect_false(isTRUE(all.equal(a, c)))
})

test_that("the null tolerates the missing values a real matrix carries", {
  d <- planted()
  d$mat[sample.int(length(d$mat), 500)] <- NA
  expect_no_error(blood_cor_null_quantile(d$mat, d$index, n_perm = 20))
})

test_that("permuting the index destroys a planted association", {
  d <- planted()
  cut <- blood_cor_null_quantile(d$mat, d$index, n_perm = 60)
  # Protein 1 correlates at ~1 with the real index; the null never sees that.
  expect_lt(cut, 0.9)
})
