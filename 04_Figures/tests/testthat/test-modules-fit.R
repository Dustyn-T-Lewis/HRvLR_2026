pacman::p_load(here, testthat, WGCNA)
source(here("04_Figures", "05_test", "modules", "a_script", "functions", "fit.R"))

test_that("fit_modules is deterministic and drops grey", {
  set.seed(1)
  n_s <- 40
  block <- 50
  k <- 4
  latent <- matrix(rnorm(n_s * k), n_s, k)
  mat <- do.call(cbind, lapply(seq_len(k), function(j) {
    latent[, j] %o% rep(1, block) + matrix(rnorm(n_s * block, sd = 0.3), n_s, block)
  }))
  colnames(mat) <- paste0("p", seq_len(ncol(mat)))
  rownames(mat) <- paste0("s", seq_len(n_s))
  meta <- data.frame(sample_id = rownames(mat))

  a <- fit_modules(t(mat), meta)
  b <- fit_modules(t(mat), meta)
  expect_identical(a$colors, b$colors)
  expect_false("MEgrey" %in% colnames(a$eigengenes))
  expect_equal(nrow(a$eigengenes), 40)
  expect_gte(ncol(a$eigengenes), 1)
})
