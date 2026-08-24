# Map a proteome change onto a phenotype, continuously.
#
# Each subject contributes one column: how much a feature moved over a window.
# Regressing that on the subject's adaptation asks whether the two track each
# other, without cutting anyone into a group. Dichotomising would throw away
# about a third of the effective sample (Cohen 1983) and invite the circularity
# a median split on a composite creates, since the split then separates its own
# ingredients whether or not the composite means anything.
#
# One row per subject means no repeated measures inside the fit, so no blocking
# and no duplicateCorrelation: the within-subject structure is spent forming
# the difference. That is also why Baseline has no place here. A baseline
# association compares levels between people, which is a different question
# from whether a change tracks a change.

pacman::p_load(here, dplyr, tidyr, tibble, readr, limma)

source(here("functions", "feature_levels.R"))

WINDOWS <- c(training = "T2 - T1", acute = "T3 - T2")

window_timepoints <- function(window) {
  switch(window,
    training = c("T1", "T2"),
    acute = c("T2", "T3"),
    stop("unknown window: ", window)
  )
}

sample_metadata <- function() {
  dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
  m <- as.data.frame(dal$metadata)
  tibble(
    sample_id = m$Col_ID,
    subject = m$Subject_ID,
    timepoint = m$Timepoint
  )
}

# Features by subject, holding the later timepoint minus the earlier one. A
# subject missing either timepoint is dropped rather than half-counted, which
# is what removes S28 from the training window and S29 from neither.
subject_change <- function(mat, window, meta = sample_metadata()) {
  tp <- window_timepoints(window)
  keyed <- meta |>
    filter(.data$timepoint %in% tp, .data$sample_id %in% colnames(mat))
  complete <- keyed |>
    count(.data$subject) |>
    filter(.data$n == 2L) |>
    pull(.data$subject)
  if (!length(complete)) {
    stop("no subject has both timepoints for window: ", window)
  }
  early <- keyed |>
    filter(.data$timepoint == tp[1], .data$subject %in% complete)
  late <- keyed |>
    filter(.data$timepoint == tp[2], .data$subject %in% complete)
  early <- early[match(complete, early$subject), ]
  late <- late[match(complete, late$subject), ]
  out <- mat[, late$sample_id, drop = FALSE] -
    mat[, early$sample_id, drop = FALSE]
  colnames(out) <- complete
  out
}

# limma across all features at once against one continuous predictor. The
# moderated variance is the reason to use it over a per-feature lm at n = 14:
# it borrows strength across features instead of trusting each one's own
# residual.
associate <- function(feat, y) {
  shared <- intersect(colnames(feat), names(y))
  y <- y[shared]
  keep <- !is.na(y)
  y <- y[keep]
  feat <- feat[, shared[keep], drop = FALSE]
  design <- stats::model.matrix(~y)
  fit <- limma::eBayes(limma::lmFit(feat, design))
  res <- limma::topTable(
    fit,
    coef = 2, number = Inf, adjust.method = "BH", sort.by = "none"
  ) |>
    tibble::rownames_to_column("feature") |>
    transmute(
      feature = .data$feature, slope = .data$logFC, t = .data$t,
      p = .data$P.Value, bh = .data$adj.P.Val
    )
  attr(res, "n") <- length(y)
  res
}

phenotype_table <- function() {
  read_csv(here("00_input", "c_data", "phenotype.csv"), show_col_types = FALSE)
}

phenotype_vector <- function(pheno, name) {
  setNames(pheno[[name]], pheno$subject)
}

# The three levels the project already fits, all keyed by the 45 sample ids so
# subject_change() treats them identically.
feature_matrices <- function() {
  list(
    proteins = protein_matrix(),
    modules = module_matrix(),
    pathways = pathway_matrix()
  )
}
