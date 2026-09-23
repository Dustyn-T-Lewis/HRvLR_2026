# Map the proteome onto a phenotype, continuously.
#
# Each subject contributes one column: a feature's value at a timepoint, or how
# much it moved over a window. Regressing that on the subject's adaptation asks
# whether the two track each other without cutting anyone into a group.
# Dichotomising would throw away about a third of the effective sample
# (Cohen 1983) and invite the circularity a median split on a composite
# creates, since the split then separates its own ingredients whether or not
# the composite means anything.
#
# One row per subject means no repeated measures inside the fit, so no blocking
# and no duplicateCorrelation: for a change window the within-subject structure
# is spent forming the difference, and for a level window only one timepoint
# enters.

pacman::p_load(here, dplyr, tidyr, tibble, readr, limma)

# Six views of the same proteome. The three levels ask a between-person
# question: do people whose module sits higher at this timepoint adapt more.
# The three changes ask a within-person one: does a shift track a shift. They
# are different questions and are reported apart, never pooled.
#
# The three levels are close to one question rather than three, because a
# subject's proteome at T1, T2 and T3 is largely the same proteome.
WINDOWS <- c(
  T1 = "level at T1",
  T2 = "level at T2",
  T3 = "level at T3",
  training = "T2 - T1",
  acute = "T3 - T2",
  total = "T3 - T1"
)

LEVEL_WINDOWS <- c("T1", "T2", "T3")

window_timepoints <- function(window) {
  switch(window,
    T1 = "T1",
    T2 = "T2",
    T3 = "T3",
    training = c("T1", "T2"),
    acute = c("T2", "T3"),
    total = c("T1", "T3"),
    stop("unknown window: ", window)
  )
}

sample_metadata <- function() {
  dal <- readRDS(here(
    "01_Preprocess", "02_Normalization", "c_data", "DAList_normalized.rds"
  ))
  m <- as.data.frame(dal$metadata)
  tibble(
    sample_id = m$Col_ID,
    subject = m$Subject_ID,
    arm = m$Group,
    timepoint = m$Timepoint
  )
}

# Features by subject for one window: the value at a timepoint, or the later
# timepoint minus the earlier one. A subject missing any timepoint the window
# needs is dropped rather than half-counted, which is what removes S28 from
# every window touching T2 or T3 and S29 from every window touching T1.
subject_window <- function(mat, window, meta = sample_metadata()) {
  tp <- window_timepoints(window)
  keyed <- meta |>
    filter(.data$timepoint %in% tp, .data$sample_id %in% colnames(mat))
  complete <- keyed |>
    count(.data$subject) |>
    filter(.data$n == length(tp)) |>
    pull(.data$subject)
  if (!length(complete)) {
    stop("no subject has every timepoint for window: ", window)
  }
  at <- function(t) {
    rows <- keyed |> filter(.data$timepoint == t, .data$subject %in% complete)
    rows <- rows[match(complete, rows$subject), ]
    mat[, rows$sample_id, drop = FALSE]
  }
  out <- if (length(tp) == 1L) at(tp) else at(tp[2]) - at(tp[1])
  colnames(out) <- complete
  out
}

# limma across all features at once against one continuous predictor. The
# moderated variance is the reason to use it over a per-feature lm at n = 14:
# it borrows strength across features instead of trusting each one's own
# residual.
#
# `adjust` is a subject-by-covariate matrix entering the design beside the
# phenotype, leaving coef 2 as the phenotype slope. Biopsy composition is the
# intended caller: at T2 the myofibre fraction correlates -0.81 with change in
# whole-muscle CSA and the blood fraction +0.71, so a level-window association
# with that phenotype is partly a statement about what the needle collected.
associate <- function(feat, y, adjust = NULL) {
  shared <- intersect(colnames(feat), names(y))
  y <- y[shared]
  keep <- !is.na(y)
  y <- y[keep]
  feat <- feat[, shared[keep], drop = FALSE]
  design <- if (is.null(adjust)) {
    stats::model.matrix(~y)
  } else {
    cov <- adjust[colnames(feat), , drop = FALSE]
    stats::model.matrix(~ y + cov)
  }
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
