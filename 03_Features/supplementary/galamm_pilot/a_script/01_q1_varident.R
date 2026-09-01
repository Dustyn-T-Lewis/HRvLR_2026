#!/usr/bin/env Rscript
# Q1: does a per-timepoint residual variance change any contrast?
#
# shared_hlm.R fits one sigma across all three timepoints, but T3 biopsies
# carry roughly twice the blood of T1 and T2 and the rise differs by arm
# (arm x T3 b = -1.21, p = 0.032), so a single sigma pools an inflated
# timepoint into the error term of every contrast. Each of the 931
# complete-case proteins is fitted twice with nlme::lme, homoscedastic and
# varIdent(~ 1 | timepoint), same fixed effects
# (group * timepoint + blood_index), random intercept per subject, REML,
# so the 2-df likelihood ratio tests the variance structure alone. nlme
# answers this cheaper than galamm and identically (PREREG.md).
#
# Read 01_q1_variance.csv for per-protein sigma ratios and the LRT, and
# 01_q1_contrasts.csv for the six shared_hlm contrasts under both fits,
# BH within model x contrast. 01_q1_lmertest.csv is the existing
# lmerTest spec on the same proteins; it uses Satterthwaite df where lme
# uses containment, so part of any p difference is the df approximation.
# A null looks like both median ratios inside [0.8, 1.25] and no contrast
# gaining a BH q < 0.05 protein the homoscedastic fit lacked.

pacman::p_load(here, dplyr, tidyr, purrr, tibble, readr, nlme, emmeans)
source(here("functions", "shared_hlm.R"))
source(here("functions", "blood_index_model.R"))
source(here(
  "03_Features", "supplementary", "galamm_pilot", "a_script",
  "pilot_helpers.R"
))

set.seed(42)
OUT <- here("03_Features", "supplementary", "galamm_pilot", "c_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

inp <- pilot_data()
mat <- inp$mat
meta <- inp$meta
CORES <- max(1L, parallel::detectCores() - 2L)

lme_contrasts <- function(fit, d, label) {
  emm <- emmeans(fit, ~ group * timepoint, data = d)
  con <- as.data.frame(summary(contrast(emm, hlm_contrast_weights(emm))))
  tibble(
    model = label, contrast = con$contrast, estimate = con$estimate,
    se = con$SE, df = con$df, t = con$t.ratio, p = con$p.value
  )
}

fit_one <- function(i) {
  d <- data.frame(y = mat[i, ], meta)
  m0 <- tryCatch(
    lme(y ~ group * timepoint + blood_index,
      random = ~ 1 | subject, data = d, method = "REML"
    ),
    error = function(e) NULL
  )
  m1 <- if (!is.null(m0)) {
    tryCatch(
      update(m0, weights = varIdent(form = ~ 1 | timepoint)),
      error = function(e) NULL
    )
  }
  feature <- rownames(mat)[i]
  if (is.null(m0) || is.null(m1)) {
    return(list(variance = tibble(
      feature = feature, ratio_t2 = NA_real_, ratio_t3 = NA_real_,
      lrt_p = NA_real_,
      status = if (is.null(m0)) "homoscedastic_failed" else "varident_failed"
    )))
  }
  ratios <- varident_ratios(m1)
  list(
    variance = tibble(
      feature = feature, ratio_t2 = ratios["T2"], ratio_t3 = ratios["T3"],
      lrt_p = anova(m0, m1)$`p-value`[2], status = "ok"
    ),
    contrasts = bind_rows(
      lme_contrasts(m0, d, "homoscedastic"),
      lme_contrasts(m1, d, "varident")
    ) |> mutate(feature = feature, .before = 1)
  )
}

fits <- parallel::mclapply(seq_len(nrow(mat)), fit_one, mc.cores = CORES)

variance <- bind_rows(map(fits, "variance"))
contrasts_tbl <- bind_rows(map(fits, "contrasts")) |>
  mutate(bh = p.adjust(p, "BH"), .by = c("model", "contrast"))

lmertest <- associate_global_hlm(
  mat, meta[, c("sample", "subject", "group", "timepoint")],
  cores = CORES
)

write_csv(variance, file.path(OUT, "01_q1_variance.csv"))
write_csv(contrasts_tbl, file.path(OUT, "01_q1_contrasts.csv"))
write_csv(lmertest, file.path(OUT, "01_q1_lmertest.csv"))

ok <- variance$status == "ok"
hits <- contrasts_tbl |>
  summarise(bh05 = sum(bh < 0.05, na.rm = TRUE), .by = c("model", "contrast"))
message(sprintf(
  "converged %d/%d | median ratios T2/T1 %.3f, T3/T1 %.3f | LRT p<.05: %d",
  sum(ok), nrow(mat),
  median(variance$ratio_t2[ok]), median(variance$ratio_t3[ok]),
  sum(variance$lrt_p[ok] < 0.05)
))
message(paste(capture.output(print(as.data.frame(
  pivot_wider(hits, names_from = "model", values_from = "bh05")
))), collapse = "\n"))
