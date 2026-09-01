source(here::here(
  "03_Features", "supplementary", "mcsa_axis", "a_script",
  "mcsa_helpers.R"
))

fake_pheno <- function() {
  data.frame(
    subject = paste0(rep(c("HR_S", "LR_S"), each = 4), 1:8),
    group_arm = rep(c("HR", "LR"), each = 4),
    comp_hypertrophy = c(10, 8, 6, 4, 2, 0, -2, -4),
    d_fcsa_I = c(3, 2, 1, 0, -1, -2, -3, -4),
    d_fcsa_II = c(6, 4, 2, 0, -2, -4, -6, -8),
    d_mcsa = c(1, -1, 2, -2, 1, -1, 2, -2)
  )
}

fake_design <- function(n = 10, effect = 1) {
  mcsa <- seq(-2, 2, length.out = n)
  list(
    y = rbind(
      signal = effect * mcsa + stats::rnorm(n, sd = 0.2),
      noise = stats::rnorm(n)
    ),
    blood = stats::rnorm(n),
    mcsa = mcsa,
    subject = paste0(rep(c("HR_S", "LR_S"), length.out = n), seq_len(n))
  )
}

test_that("fibre_axis averages the two fibre measures on a common scale", {
  pheno <- fake_pheno()
  axis <- fibre_axis(pheno)
  expect_equal(mean(axis), 0, tolerance = 1e-10)
  expect_equal(stats::sd(axis), 1, tolerance = 1e-10)
  expect_equal(
    stats::cor(axis, pheno$d_fcsa_I, method = "spearman"), 1,
    tolerance = 1e-10
  )
})

test_that("fibre_axis ignores the scale each measure came in on", {
  pheno <- fake_pheno()
  rescaled <- transform(pheno, d_fcsa_II = d_fcsa_II * 1000)
  expect_equal(fibre_axis(pheno), fibre_axis(rescaled), tolerance = 1e-10)
})

test_that("composite_weights recovers a composite it was built from", {
  pheno <- fake_pheno()
  z <- scale(pheno[, CSA_ITEMS])
  pheno$comp_hypertrophy <- 3 * z[, "d_fcsa_I"] + 1 * z[, "d_mcsa"] +
    withr::with_seed(4, stats::rnorm(8, sd = 0.01))

  weights <- composite_weights(pheno)
  expect_equal(weights$item, CSA_ITEMS)
  expect_equal(weights$weight[weights$item == "d_mcsa"], 1, tolerance = 0.02)
  expect_gt(weights$r2_full[1], 0.999)
  expect_lt(weights$r2_without_mcsa[1], weights$r2_full[1])
})

test_that("arm_auc is 1 for perfect separation and 0.5 for none", {
  arm <- rep(c("HR", "LR"), each = 4)
  expect_equal(arm_auc(c(4, 3, 2, 1, 0, -1, -2, -3), arm), 1)
  expect_equal(arm_auc(rep(1, 8), arm), 0.5)
  expect_equal(arm_auc(c(1, 1, 1, 1, 2, 2, 2, 2), arm), 0)
})

test_that("arm_auc drops a subject whose value is missing", {
  arm <- rep(c("HR", "LR"), each = 4)
  expect_equal(arm_auc(c(4, 3, 2, NA, 0, -1, -2, -3), arm), 1)
})

test_that("config_design differences both response and covariate for delta", {
  mat <- matrix(
    c(1, 2, 4, 8, 3, 5),
    nrow = 1, dimnames = list("p", c("a1", "b1", "a2", "b2", "a3", "b3"))
  )
  meta <- data.frame(
    sample = colnames(mat), subject = c("HR_S1", "LR_S2"),
    timepoint = rep(c("T1", "T2", "T3"), each = 2),
    blood_index = c(1, 2, 5, 8, 0, 0)
  )
  pheno <- data.frame(subject = c("HR_S1", "LR_S2"), d_mcsa = c(0.5, -0.5))

  delta <- config_design(mat, meta, pheno, "delta")
  expect_equal(delta$subject, c("HR_S1", "LR_S2"))
  expect_equal(as.numeric(delta$y), c(3, 6))
  expect_equal(delta$blood, c(4, 6))
  expect_equal(delta$mcsa, c(0.5, -0.5))
})

test_that("config_design drops a subject and keeps the phenotype aligned", {
  mat <- matrix(1:4, nrow = 1, dimnames = list("p", c("a1", "b1", "a2", "b2")))
  meta <- data.frame(
    sample = colnames(mat), subject = c("HR_S1", "LR_S2"),
    timepoint = rep(c("T1", "T2"), each = 2), blood_index = 1:4
  )
  pheno <- data.frame(subject = c("HR_S1", "LR_S2"), d_mcsa = c(0.5, -0.5))

  kept <- config_design(mat, meta, pheno, "T1", drop = "LR_S2")
  expect_equal(kept$subject, "HR_S1")
  expect_equal(kept$mcsa, 0.5)
  expect_equal(ncol(kept$y), 1L)
})

test_that("mcsa_scan separates a planted slope from noise", {
  withr::with_seed(1, {
    scan <- mcsa_scan(fake_design(n = 12, effect = 1))
  })
  expect_equal(scan$feature, c("signal", "noise"))
  expect_gt(scan$slope[1], 0.7)
  expect_lt(scan$p[1], 1e-6)
  expect_gt(scan$p[2], 0.05)
  expect_true(all(scan$bh >= scan$p))
})

test_that("mcsa_scan reports no slope when the outcome carries none", {
  withr::with_seed(1, {
    scan <- mcsa_scan(fake_design(n = 12, effect = 0))
  })
  expect_gt(min(scan$p), 0.05)
  expect_gt(min(scan$bh), 0.05)
})

test_that("mcsa_permutation floors p for the planted slope, not for noise", {
  design <- withr::with_seed(1, fake_design(n = 12, effect = 1))
  observed <- mcsa_scan(design)
  perm <- mcsa_permutation(design, observed, n_perm = 50, seed = 3)

  expect_equal(perm$feature, c("signal", "noise"))
  expect_equal(perm$emp_p[1], 1 / 51, tolerance = 1e-8)
  expect_equal(perm$n_ge[1], 0L)
  expect_gt(perm$emp_p[2], 1 / 51)
  expect_true(all(perm$emp_p <= 1))
})

test_that("mcsa_permutation handles a single surviving feature", {
  design <- withr::with_seed(1, fake_design(n = 12, effect = 1))
  observed <- mcsa_scan(design)[1, ]
  perm <- mcsa_permutation(design, observed, n_perm = 20, seed = 3)

  expect_equal(nrow(perm), 1L)
  expect_equal(perm$emp_p, 1 / 21, tolerance = 1e-8)
})

test_that("mcsa_permutation is reproducible from its seed", {
  design <- withr::with_seed(1, fake_design(n = 12, effect = 0.4))
  observed <- mcsa_scan(design)
  expect_equal(
    mcsa_permutation(design, observed, n_perm = 20, seed = 7),
    mcsa_permutation(design, observed, n_perm = 20, seed = 7)
  )
})

test_that("survivor_checks tells a within-arm slope from an arm difference", {
  design <- withr::with_seed(1, fake_design(n = 12, effect = 1))
  within_arm <- survivor_checks(design, "signal")
  expect_gt(within_arm$loo_t_min, 3)
  expect_gt(within_arm$rho_hr, 0.8)
  expect_gt(within_arm$rho_lr, 0.8)
  expect_gt(within_arm$p_arm_term, 0.05)

  arm_only <- design
  arm_only$y["signal", ] <- as.numeric(startsWith(design$subject, "HR")) +
    withr::with_seed(2, stats::rnorm(12, sd = 0.05))
  expect_lt(survivor_checks(arm_only, "signal")$p_arm_term, 0.01)
})
