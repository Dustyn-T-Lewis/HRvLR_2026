# Differential variability: does a protein's spread across subjects differ
# between HR and LR, not just its mean? DEP (01_run_dep.R) tests the mean
# shift; this is a second, independent lens on the same nine contrasts,
# looking for a discriminating signal a mean-shift test cannot see.
#
# Phipson & Oshlack 2014, Genome Biology 15:465 -- DiffVar, inspired by
# Levene's test: fit a linear model to each protein's absolute deviation
# from its group median, then moderate with the same limma empirical-Bayes
# engine DEP already uses. Built for bulk methylation arrays, not single-cell
# data, and reuses the exact zero-intercept design and contrast matrix
# feature_design() already builds for DEP -- no new design-building code.
#
# Two real limitations, not shortcuts:
#   - missMethyl::varFit() has no block or correlation argument (confirmed
#     against the package source and Bioconductor support threads), so
#     unlike every other test in this pipeline, this one does NOT carry the
#     subject-blocking that duplicateCorrelation supplies elsewhere. A hit
#     here is not yet adjusted for repeated measures the way a DEP or fry
#     hit is; treat it as a lead to look at, not a confirmed result.
#   - varFit()'s range check (min/max of the input) errors on any NA, so it
#     cannot take the primary non-imputed matrix DEP uses. This runs on the
#     missForest arm instead, the same complete-data requirement
#     02_Pathways/01_run_pathways.R's fry call already carries -- an
#     asymmetry every reader of this table has to state alongside it.

pacman::p_load(here, dplyr, tibble, openxlsx)
source(here("functions", "feature_contrasts.R"))
source(here("functions", "feature_levels.R"))
source(here("functions", "pred_features.R"))
source(here(
  "03_Features", "01_Proteins", "a_script", "pi_permutation.R"
))

mat <- as.matrix(readRDS(pred_paths()$dalist)$data)
meta <- feature_metadata()
parts <- feature_design(mat, meta)

vfit <- missMethyl::varFit(
  mat,
  design = parts$design, coef = seq_len(ncol(parts$design))
)
vfit2 <- missMethyl::contrasts.varFit(vfit, contrasts = parts$contrasts)

results <- bind_rows(lapply(colnames(parts$contrasts), function(ct) {
  missMethyl::topVar(vfit2,
    coef = ct, number = nrow(mat), sort = FALSE
  ) |>
    tibble::rownames_to_column("feature") |>
    transmute(
      contrast = ct, feature = .data$feature,
      LogVarRatio = .data$LogVarRatio, t = .data$t,
      p = .data$P.Value, bh = .data$Adj.P.Value
    )
}))

write.xlsx(
  list(diffvar = results),
  here("03_Features", "01_Proteins", "c_data", "08_diffvar.xlsx")
)

summary_tbl <- results |>
  summarise(
    n_tested = dplyr::n(), n_nominal = sum(.data$p < 0.05),
    n_bh = sum(.data$bh < 0.05), min_bh = min(.data$bh),
    .by = "contrast"
  )
print(as.data.frame(summary_tbl), digits = 3)

hits <- filter(results, .data$bh < 0.05)

if (nrow(hits) == 0) {
  cat("\nNo BH<0.05 differential-variability hits in any contrast.\n")
  cat("Permutation check skipped: nothing to confirm.\n")
} else {
  cat(sprintf(
    "\n%d differential-variability hit(s) at BH<0.05 -- confirming against a",
    nrow(hits)
  ))
  cat(" subject-label permutation null (200 shuffles), matching the\n")
  cat("04_pi_permutation.R pattern.\n")
  print(as.data.frame(hits), digits = 3)

  N_PERM <- 200
  observed_min_bh <- min(results$bh)

  fit_permuted_min_bh <- function(perm_meta) {
    design <- stats::model.matrix(~ 0 + group, perm_meta)
    colnames(design) <- levels(perm_meta$group)
    cm <- limma::makeContrasts(contrasts = HRVLR_CONTRASTS, levels = design)
    colnames(cm) <- trimws(sub("=.*$", "", HRVLR_CONTRASTS))
    pvfit <- missMethyl::varFit(
      mat,
      design = design, coef = seq_len(ncol(design))
    )
    pvfit2 <- missMethyl::contrasts.varFit(pvfit, contrasts = cm)
    min(vapply(colnames(cm), function(ct) {
      min(missMethyl::topVar(
        pvfit2,
        coef = ct, number = nrow(mat), sort = FALSE
      )$Adj.P.Value)
    }, numeric(1)))
  }

  set.seed(42)
  meta_ord <- meta[match(colnames(mat), meta$sample_id), ]
  null_min_bh <- vapply(seq_len(N_PERM), function(i) {
    suppressWarnings(fit_permuted_min_bh(permute_arm_labels(meta_ord)))
  }, numeric(1))

  emp_p <- (sum(null_min_bh <= observed_min_bh) + 1) / (N_PERM + 1)
  cat(sprintf(
    "\nObserved smallest BH across all contrasts: %.4f\n", observed_min_bh
  ))
  cat(sprintf(
    "Permuted median smallest BH (n=%d): %.4f\n", N_PERM, median(null_min_bh)
  ))
  cat(sprintf("Empirical p (permutation as extreme or more): %.4f\n", emp_p))

  write.xlsx(
    list(
      diffvar = results, hits = hits,
      permutation = tibble(
        n_perm = N_PERM, observed_min_bh = observed_min_bh,
        null_median_min_bh = median(null_min_bh), emp_p = emp_p
      )
    ),
    here(
      "03_Features", "01_Proteins", "c_data", "08_diffvar.xlsx"
    )
  )
}
