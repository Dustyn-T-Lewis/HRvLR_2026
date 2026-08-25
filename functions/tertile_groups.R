# Three responder groups, and the one question they can answer that a
# two-group split and a linear fit cannot.
#
# Cutting 16 subjects into thirds buys no power. Simulated on this cohort's own
# composite, a 5-v-5 extremes contrast detects a feature linear in the trait
# less often than the 8-v-8 median split (30.9% vs 35.9% at beta = 0.5), and
# both lose to regressing on all 16 (43.7%). The gain in separation, 1.78 to
# 2.23 SD, is cancelled almost exactly by the six subjects it costs.
#
# What three groups add is shape. A linear fit is blind by construction to a
# non-monotonic relationship: if mid-responders sit above or below both
# extremes, every model in this project so far misses it. Splitting the
# three-group effect into an ordered linear trend and a quadratic deviation
# separates "tracks the trait" from "the middle is different", and only the
# second is new.
#
# Groups are cut on a phenotype, never on the proteome, so nothing here is
# circular. The middle group is retained and modelled rather than discarded:
# dropping it would throw away the only subjects that carry the shape.

# matrixStats is called qualified, never attached: matrixStats::count() masks
# dplyr::count(), which subject_window() relies on.
pacman::p_load(here, dplyr, tidyr, tibble, readr, limma)

source(here("functions", "feature_levels.R"))
source(here("functions", "association.R"))

GROUPS <- c("LR", "MR", "HR")

# Equal thirds where n divides by three, otherwise the extra subjects go to the
# middle so the extremes stay clean. At n = 16 that gives 5 / 6 / 5.
tertile_split <- function(x, subject) {
  keep <- !is.na(x)
  x <- x[keep]
  subject <- subject[keep]
  n <- length(x)
  n_edge <- floor(n / 3)
  r <- rank(x, ties.method = "first")
  group <- dplyr::case_when(
    r <= n_edge ~ "LR",
    r > n - n_edge ~ "HR",
    TRUE ~ "MR"
  )
  tibble(subject = subject, group = factor(group, levels = GROUPS))
}

# Orthogonal polynomial contrasts on the ordered factor: coefficient 2 is the
# linear trend across LR < MR < HR, coefficient 3 the quadratic deviation that
# a two-group or linear model cannot express. Reporting both from one fit keeps
# them on the same residual variance.
tertile_design <- function(group) {
  g <- factor(as.character(group), levels = GROUPS, ordered = TRUE)
  stats::contrasts(g) <- stats::contr.poly(3)
  design <- stats::model.matrix(~g)
  colnames(design) <- c("intercept", "linear", "quadratic")
  design
}

# One fit, three read-outs: the ordered trend, the non-monotonic deviation, and
# the omnibus F for any group difference at all. BH is applied within each
# term, never pooled across the three.
#
# `min_per_group` is not optional tidying. The protein matrix is the
# non-imputed one and is 12% missing, lmFit drops NAs row by row, and a change
# window compounds it by needing both timepoints. In the acute window 141
# proteins carry fewer than five observations. Without this filter a feature
# measured in one subject of one group still returns a group contrast, and the
# first quadratic hit this sweep produced was exactly that: latexin, with a
# single observation behind its high group.
fit_tertile <- function(feat, group_tbl, min_per_group = 3L) {
  shared <- intersect(colnames(feat), group_tbl$subject)
  g <- group_tbl[match(shared, group_tbl$subject), ]
  feat <- feat[, shared, drop = FALSE]

  observed <- vapply(GROUPS, function(grp) {
    cols <- which(g$group == grp)
    rowSums(!is.na(feat[, cols, drop = FALSE]))
  }, numeric(nrow(feat)))
  keep <- matrixStats::rowMins(observed) >= min_per_group
  dropped <- sum(!keep)
  feat <- feat[keep, , drop = FALSE]
  if (!nrow(feat)) {
    stop("no feature has ", min_per_group, " observations in every group")
  }

  design <- tertile_design(g$group)
  fit <- limma::eBayes(limma::lmFit(feat, design))

  term <- function(coef, label) {
    limma::topTable(
      fit,
      coef = coef, number = Inf, adjust.method = "BH", sort.by = "none"
    ) |>
      tibble::rownames_to_column("feature") |>
      transmute(
        term = label, feature = .data$feature, estimate = .data$logFC,
        t = .data$t, p = .data$P.Value, bh = .data$adj.P.Val
      )
  }
  omnibus <- limma::topTable(
    fit,
    coef = 2:3, number = Inf, adjust.method = "BH", sort.by = "none"
  ) |>
    tibble::rownames_to_column("feature") |>
    transmute(
      term = "omnibus", feature = .data$feature,
      estimate = NA_real_, t = .data$F,
      p = .data$P.Value, bh = .data$adj.P.Val
    )

  res <- bind_rows(term(2, "linear"), term(3, "quadratic"), omnibus)
  attr(res, "n") <- ncol(feat)
  attr(res, "sizes") <- paste(table(g$group), collapse = "/")
  attr(res, "n_features") <- nrow(feat)
  attr(res, "n_dropped") <- dropped
  res
}
