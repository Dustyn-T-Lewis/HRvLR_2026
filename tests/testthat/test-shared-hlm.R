test_that("fit_hlm reproduces a direct lmerTest REML fit", {
  source(here::here("functions", "shared_hlm.R"))
  pacman::p_load(lmerTest, lme4)

  set.seed(1)
  d <- expand.grid(subject = factor(1:12), Timepoint = factor(c("T1", "T2")))
  d$Group <- factor(rep(c("HR", "LR"), each = 12))
  d$value <- rnorm(nrow(d)) + as.integer(d$Timepoint) * 0.5

  got <- fit_hlm(d, response = "value")
  ref <- lmerTest::lmer(
    value ~ Group * Timepoint + (1 | subject),
    data = d, REML = TRUE
  )
  expect_s4_class(got, "lmerModLmerTest")
  expect_equal(lme4::fixef(got), lme4::fixef(ref), tolerance = 1e-8)
  # Satterthwaite p-values survive (this is why lmerTest, not plain lme4)
  expect_true("Pr(>|t|)" %in% colnames(summary(got)$coefficients))
})

.hlm_toy <- function(seed = 1) {
  set.seed(seed)
  subj <- sprintf("S%02d", 1:16)
  grp <- rep(c("HR", "LR"), each = 8)
  meta <- expand.grid(
    subject = subj, timepoint = c("T1", "T2", "T3"),
    stringsAsFactors = FALSE
  )
  meta$group <- grp[match(meta$subject, subj)]
  meta$sample <- paste(meta$subject, meta$timepoint, sep = "_")
  meta$group <- factor(meta$group, levels = c("LR", "HR"))
  meta$timepoint <- factor(meta$timepoint, levels = c("T1", "T2", "T3"))
  base <- rnorm(16, sd = 1.5)[match(meta$subject, subj)]
  signal <- 2 * (meta$group == "HR") * (meta$timepoint == "T2")
  mat <- rbind(
    hit = base + signal + rnorm(nrow(meta), sd = 0.3),
    noise = base + rnorm(nrow(meta), sd = 0.3)
  )
  colnames(mat) <- meta$sample
  list(mat = mat, meta = meta)
}

test_that("associate_global_hlm returns the six contrasts with a T2 hit", {
  source(here::here("functions", "shared_hlm.R"))
  toy <- .hlm_toy()

  res <- suppressWarnings(associate_global_hlm(toy$mat, toy$meta))
  expect_setequal(unique(res$contrast), HLM_CONTRASTS)
  expect_named(res, c(
    "feature", "contrast", "estimate", "se", "df", "t", "p", "bh"
  ))

  # the planted HR-at-T2 bump shows up as a positive training contrast on `hit`
  tr <- res[res$feature == "hit" & res$contrast == "training", ]
  expect_gt(tr$estimate, 0)
  expect_lt(tr$p, 0.05)
  # and the noise protein's per-timepoint contrasts are not significant
  t2_noise <- res[res$feature == "noise" & res$contrast == "T2", ]
  expect_gt(t2_noise$p, 0.05)
})

test_that("associate_global_hlm returns NA rows for an unfittable protein", {
  source(here::here("functions", "shared_hlm.R"))
  toy <- .hlm_toy()
  toy$mat["noise", ] <- NA_real_

  res <- suppressWarnings(associate_global_hlm(toy$mat, toy$meta))
  na_rows <- res[res$feature == "noise", ]
  expect_equal(nrow(na_rows), length(HLM_CONTRASTS))
  expect_true(all(is.na(na_rows$p)))
})

test_that("varpart_global fractions sum to one and rank subject over group", {
  source(here::here("functions", "shared_hlm.R"))
  toy <- .hlm_toy()

  vp <- suppressWarnings(varpart_global(toy$mat, toy$meta))
  frac <- vp[, c("subject", "group", "timepoint", "Residuals")]
  expect_equal(rowSums(frac), c(1, 1), tolerance = 1e-6, ignore_attr = TRUE)
  expect_gt(vp$subject[vp$feature == "noise"], vp$group[vp$feature == "noise"])
})
