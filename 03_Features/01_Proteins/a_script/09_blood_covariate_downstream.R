#!/usr/bin/env Rscript
# Does adjusting for blood change anything at module, pathway or classification
# level?
#
# 08 put the covariate in the protein model. The feature levels sit on top of
# that model and the classification screen sits on the same matrix, so both
# should be asked the same question rather than assumed unaffected.
#
# Modules and pathways go through fit_feature_contrasts(adjust=), which enters
# blood as one design column carrying zero weight in every contrast, so the nine
# estimands are unchanged and only the residual moves.
#
# Classification residualises each protein on the blood index first. That is not
# outcome leakage: the index is built from five haemoglobins the contaminant
# filter removed, and it never sees the arm label. The permutation null is
# recomputed on the residualised matrix so any variance the adjustment removes
# is removed from observed and null alike.

pacman::p_load(here, dplyr, tibble, purrr, readr, mixOmics, pROC, openxlsx)
source(here("functions", "feature_contrasts.R"))
source(here("functions", "feature_levels.R"))
source(here("functions", "blood_index_model.R"))
source(here("functions", "shared_utils.R"))

set.seed(42)
B <- 200
KEEP <- 50
OUT <- here(
  "03_Features", "01_Proteins", "c_data", "blood_covariate"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

blood <- blood_index_data()
covar <- stats::setNames(blood$blood_index, blood$Col_ID)

count_hits <- function(mat, label) {
  a <- fit_feature_contrasts(mat)
  b <- fit_feature_contrasts(mat, adjust = covar)
  inner_join(
    a |> summarise(
      nom = sum(p < 0.05), n_bh = sum(bh < 0.05),
      min_bh = min(bh), .by = contrast
    ),
    b |> summarise(
      nom_adj = sum(p < 0.05), n_bh_adj = sum(bh < 0.05),
      min_bh_adj = min(bh), .by = contrast
    ),
    by = "contrast"
  ) |>
    mutate(level = label, .before = 1)
}

levels_tbl <- bind_rows(
  count_hits(module_matrix(), "modules"),
  count_hits(pathway_matrix(), "pathways")
)

cat("=== modules and pathways, primary against blood-adjusted ===\n")
levels_tbl |>
  mutate(responder = contrast %in% RESPONDER_CONTRASTS) |>
  arrange(level, min_bh_adj) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

# Classification on a blood-residualised matrix.
mat <- protein_matrix()
mat <- mat[stats::complete.cases(mat), , drop = FALSE]
b <- covar[colnames(mat)]
resid_mat <- t(apply(mat, 1, function(v) stats::residuals(stats::lm(v ~ b))))
colnames(resid_mat) <- colnames(mat)

meta <- feature_metadata()
arm <- sub("_.*$", "", as.character(meta$group))
subject <- meta$subject

loso_acc <- function(x, y, subj) {
  ok <- vapply(unique(subj), function(s) {
    tr <- subj != s
    if (length(unique(y[tr])) < 2) {
      return(NA_real_)
    }
    fit <- mixOmics::splsda(
      x[tr, , drop = FALSE], y[tr],
      ncomp = 2, keepX = rep(KEEP, 2)
    )
    mean(predict(fit, x[!tr, , drop = FALSE])$class$max.dist[, 2] == y[!tr])
  }, numeric(1))
  mean(ok, na.rm = TRUE)
}

permute_arm <- function() {
  key <- distinct(tibble(subject, arm), subject, arm)
  key$arm <- sample(key$arm)
  key$arm[match(subject, key$subject)]
}

screen <- function(x, label) {
  y <- factor(arm)
  obs <- loso_acc(t(x), y, subject)
  null <- vapply(
    seq_len(B),
    function(i) loso_acc(t(x), factor(permute_arm()), subject),
    numeric(1)
  )
  baseline <- max(table(y)) / length(y)
  cat(sprintf(
    "%-22s LOSO %.3f | null median %.3f [%.3f, %.3f] | baseline %.3f | p = %.3f\n",
    label, obs, stats::median(null), stats::quantile(null, 0.025),
    stats::quantile(null, 0.975), baseline, (1 + sum(null >= obs)) / (1 + B)
  ))
  tibble(
    matrix = label, loso = obs, null_median = stats::median(null),
    baseline = baseline, p = (1 + sum(null >= obs)) / (1 + B)
  )
}

cat("\n=== can the proteome classify HR vs LR once blood is removed? ===\n")
class_tbl <- bind_rows(
  screen(mat, "unadjusted"),
  screen(resid_mat, "blood-residualised")
)

write.xlsx(
  list(feature_levels = levels_tbl, classification = class_tbl),
  file.path(OUT, "blood_covariate_downstream.xlsx")
)
cat(sprintf("\nwrote %s\n", file.path(OUT, "blood_covariate_downstream.xlsx")))
