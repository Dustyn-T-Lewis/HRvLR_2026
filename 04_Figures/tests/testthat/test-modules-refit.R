pacman::p_load(here, testthat)
source(here("04_Figures", "05_test", "modules", "a_script", "functions", "fit.R"))
source(here("04_Figures", "05_test", "modules", "a_script", "functions", "coupling.R"))
source(here("04_Figures", "05_test", "modules", "a_script", "functions", "prediction.R"))
source(here("04_Figures", "05_test", "modules", "a_script", "functions", "honest_refit.R"))

test_that("jaccard and pc1 projection behave", {
  expect_equal(jaccard(c("a", "b", "c"), c("b", "c", "d")), 2 / 4)

  set.seed(3)
  m <- matrix(rnorm(30 * 6), nrow = 30)
  fit <- fit_pc1(m[1:24, ])
  s_in <- project_pc1(fit, m[1:24, ])
  s_out <- project_pc1(fit, m[25:30, ])
  expect_length(s_out, 6)
  expect_gt(abs(cor(s_in, rowMeans(scale(m[1:24, ])))), 0.5)
})

test_that("match_modules pairs by best jaccard", {
  full <- c(p1 = "red", p2 = "red", p3 = "blue", p4 = "blue")
  train <- c(p1 = "green", p2 = "green", p3 = "yellow", p4 = "yellow")
  mm <- match_modules(full, train)
  expect_equal(mm$train[mm$full == "red"], "green")
  expect_equal(mm$jaccard[mm$full == "red"], 1)
})

test_that("run_honest_refit returns a drop per surviving module", {
  skip_on_cran()
  set.seed(5)
  subj <- paste0("S", 1:8)
  samples <- as.vector(t(outer(subj, c("T1", "T2"), paste, sep = "_")))
  ns <- length(samples)
  block <- 45
  k <- 3
  latent <- matrix(rnorm(ns * k), ns, k)
  abund <- t(do.call(cbind, lapply(seq_len(k), function(j) {
    latent[, j] %o% rep(1, block) + matrix(rnorm(ns * block, sd = 0.3), ns, block)
  })))
  rownames(abund) <- paste0("p", seq_len(nrow(abund)))
  colnames(abund) <- samples
  meta <- data.frame(
    sample_id = samples, subject = rep(subj, each = 2),
    timepoint = rep(c("T1", "T2"), 8)
  )
  y <- setNames(rep(c(0, 1), length.out = 8), subj)

  out <- run_honest_refit(abund, meta, y, "training")
  expect_true(all(c("module", "auc_insample", "auc_loso", "drop") %in% names(out)))
  expect_gte(nrow(out), 1)
})
