test_that("associate_limma matches a direct per-subject limma fit", {
  source(here::here("04_Figures", "functions", "f04_association.R"))
  pacman::p_load(limma)

  set.seed(1)
  meta <- data.frame(
    Col_ID = paste0("s", 1:24),
    Subject_ID = rep(paste0("S", 1:8), each = 3),
    Timepoint = rep(c("T1", "T2", "T3"), times = 8)
  )
  mat <- matrix(rnorm(5 * 24),
    nrow = 5,
    dimnames = list(paste0("f", 1:5), meta$Col_ID)
  )
  pheno <- data.frame(subject = paste0("S", 1:8), d_mcsa = rnorm(8))

  got <- associate_limma(mat, meta, pheno, "d_mcsa", phases = "baseline")
  # reference: baseline = T1 columns, one per subject, limma on the trait
  t1 <- meta$Col_ID[meta$Timepoint == "T1"]
  fmat <- mat[, t1]
  colnames(fmat) <- meta$Subject_ID[meta$Timepoint == "T1"]
  y <- pheno$d_mcsa[match(colnames(fmat), pheno$subject)]
  ref <- limma::eBayes(limma::lmFit(fmat, model.matrix(~y)))
  ref_tt <- limma::topTable(ref, coef = 2, number = Inf, sort.by = "none")
  expect_equal(got$beta[match(rownames(ref_tt), got$feature)], ref_tt$logFC,
    tolerance = 1e-8
  )
  expect_true(all(c("feature", "trait", "phase", "beta", "t", "p", "bh") %in% names(got)))
})

test_that("associate_hlm returns the interaction F-test and guards sparse features", {
  source(here::here("04_Figures", "functions", "f04_association.R"))
  pacman::p_load(lmerTest)

  set.seed(1)
  meta <- data.frame(
    Col_ID = paste0("s", 1:24),
    Subject_ID = rep(paste0("S", 1:8), each = 3),
    Timepoint = rep(c("T1", "T2", "T3"), times = 8)
  )
  resp <- setNames(rnorm(8), paste0("S", 1:8))
  # feature 1: trajectory depends on response; feature 2: all NA (sparse guard)
  base <- rnorm(24)
  slope <- resp[meta$Subject_ID] * (as.integer(factor(meta$Timepoint)) - 1)
  mat <- rbind(f1 = base + slope, f2 = NA_real_)
  colnames(mat) <- meta$Col_ID
  pheno <- data.frame(subject = paste0("S", 1:8), d_mcsa = resp)

  got <- associate_hlm(mat, meta, pheno, "d_mcsa")
  expect_true(all(c("feature", "trait", "p", "bh") %in% names(got)))
  # sparse feature returns NA, does not error
  expect_true(is.na(got$p[got$feature == "f2"]))
  # a real interaction is detected as a finite F-test p for f1
  expect_true(is.finite(got$p[got$feature == "f1"]))
  # reference: the interaction p equals lmerTest anova on the same data for f1
  d <- data.frame(
    score = mat["f1", ], response = resp[meta$Subject_ID],
    timepoint = factor(meta$Timepoint), subject = meta$Subject_ID
  )
  ref <- anova(lmerTest::lmer(score ~ response * timepoint + (1 | subject), d))
  expect_equal(got$p[got$feature == "f1"],
    ref["response:timepoint", "Pr(>F)"],
    tolerance = 1e-6
  )
})
