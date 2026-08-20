# Differential variability, continuous tree: does a protein's spread across
# subjects change from Training or the Acute bout, not just its mean?
# Mirrors categorical/01_Proteins/a_script/07_diffvar.R -- same method, same
# missMethyl::varFit() -> contrasts.varFit() -> topVar() chain, same two
# stated limitations (varFit carries no block/correlation argument, so this
# test does not see the subject-blocking every other test in this pipeline
# uses; and varFit's range check errors on any NA, so this runs on the
# missForest-imputed arm, not the primary non-imputed matrix) -- fit here on
# the continuous tree's own design and two pooled contrasts instead of the
# nine.

pacman::p_load(here, dplyr, tibble, openxlsx)
source(here("functions", "feature_contrasts.R"))
source(here("functions", "feature_levels.R"))
source(here("functions", "pred_features.R"))

mat <- as.matrix(readRDS(pred_paths()$dalist)$data)
meta <- feature_metadata_pooled()
parts <- feature_design_pooled(mat, meta)

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
  here("03_Analysis", "continuous", "01_Proteins", "c_data", "04_diffvar.xlsx")
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
  cat(" within-subject timepoint-label permutation null (200 shuffles).\n")
  print(as.data.frame(hits), digits = 3)

  # Permute which of a subject's own samples counts as T1/T2/T3, keeping
  # each subject's realized timepoint set fixed (S28 stays T1-only, S29
  # stays T2/T3-only). This breaks the true sample-to-timepoint
  # correspondence while preserving the repeated-measures structure and the
  # missingness pattern the observed fit already has.
  permute_timepoint_labels <- function(meta) {
    meta |>
      mutate(
        timepoint = sample(as.character(.data$timepoint)),
        .by = "subject"
      ) |>
      mutate(timepoint = factor(.data$timepoint, levels = c("T1", "T2", "T3")))
  }

  N_PERM <- 200
  observed_min_bh <- min(results$bh)

  fit_permuted_min_bh <- function(perm_meta) {
    design <- stats::model.matrix(~ 0 + timepoint, perm_meta)
    colnames(design) <- levels(perm_meta$timepoint)
    cm <- limma::makeContrasts(contrasts = POOLED_CONTRASTS, levels = design)
    colnames(cm) <- trimws(sub("=.*$", "", POOLED_CONTRASTS))
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
    suppressWarnings(fit_permuted_min_bh(permute_timepoint_labels(meta_ord)))
  }, numeric(1))

  emp_p <- (sum(null_min_bh <= observed_min_bh) + 1) / (N_PERM + 1)
  cat(sprintf(
    "\nObserved smallest BH across both contrasts: %.4f\n", observed_min_bh
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
      "03_Analysis", "continuous", "01_Proteins", "c_data", "04_diffvar.xlsx"
    )
  )
}
