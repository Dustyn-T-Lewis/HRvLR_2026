# Does a Pi count mean anything?
#
# Pi = p^|log2FC| (Xiao Eq.2) is a transformed raw p-value and controls no error
# rate, so a count of Pi < 0.05 hits has no reference on its own. BH cannot
# supply one: BH needs a statistic uniform on [0,1] under the null, and Pi is
# not - |log2FC| is itself random, driving Pi toward 1 for the many small-effect
# proteins and toward 0 for the few large ones. The null has to be built.
#
# The permutation reassigns the HR/LR label across subjects and refits. Subject
# is the unit of randomisation because the label is a subject property;
# permuting samples would break the repeated measures and manufacture a null
# easier to beat than the design allows. Timepoint travels with the sample
# untouched, so the within-arm contrasts are permuted on the same footing as the
# group ones.
#
# The consensus correlation is estimated once on the observed data and reused.
# Re-running duplicateCorrelation per permutation costs far more than it
# changes, and the correlation is a property of the repeated-measures structure,
# which the permutation leaves intact.

pacman::p_load(here, dplyr, tibble, limma)
source(here("functions", "feature_contrasts.R"))

permute_arm_labels <- function(meta) {
  timepoint <- sub("^[^_]*_", "", as.character(meta$group))
  subjects <- unique(meta$subject)
  subject_arm <- sub("_.*$", "", as.character(
    meta$group[match(subjects, meta$subject)]
  ))
  shuffled <- setNames(sample(subject_arm), subjects)
  meta$group <- factor(
    paste(shuffled[meta$subject], timepoint, sep = "_"),
    levels = GROUP_LEVELS
  )
  meta
}

count_pi_hits <- function(res, pi_thresh) {
  res |>
    mutate(pi = pi_score(.data$p, .data$logFC)) |>
    summarise(
      n_pi = sum(.data$pi < pi_thresh, na.rm = TRUE),
      .by = "contrast"
    )
}

pi_permutation_null <- function(mat, meta = feature_metadata(), n_perm = 200,
                                seed = 42, pi_thresh = 0.05) {
  parts <- feature_design(mat, meta)
  observed <- count_pi_hits(fit_feature_contrasts(mat, meta), pi_thresh)

  fit_permuted <- function(perm_meta) {
    design <- stats::model.matrix(~ 0 + group, perm_meta)
    colnames(design) <- levels(perm_meta$group)
    cm <- limma::makeContrasts(
      contrasts = HRVLR_CONTRASTS, levels = design
    )
    colnames(cm) <- trimws(sub("=.*$", "", HRVLR_CONTRASTS))
    fit <- limma::lmFit(mat, design,
      block = parts$block, correlation = parts$correlation
    )
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
    bind_rows(lapply(colnames(cm), function(ct) {
      tt <- limma::topTable(fit2,
        coef = ct, number = Inf, sort.by = "none"
      )
      tibble(contrast = ct, logFC = tt$logFC, p = tt$P.Value)
    }))
  }

  set.seed(seed)
  meta <- meta[match(colnames(mat), meta$sample_id), ]
  null_counts <- bind_rows(lapply(seq_len(n_perm), function(i) {
    counts <- count_pi_hits(
      suppressWarnings(fit_permuted(permute_arm_labels(meta))), pi_thresh
    )
    counts$perm <- i
    counts
  }))

  null_counts |>
    summarise(
      null_median = stats::median(.data$n_pi),
      null_p95 = stats::quantile(.data$n_pi, 0.95, names = FALSE),
      n_ge = sum(.data$n_pi >= observed$n_pi[match(
        .data$contrast[1], observed$contrast
      )]),
      .by = "contrast"
    ) |>
    left_join(rename(observed, observed = "n_pi"), by = "contrast") |>
    mutate(
      n_perm = n_perm,
      emp_p = (.data$n_ge + 1) / (n_perm + 1)
    ) |>
    select(
      "contrast", "observed", "null_median", "null_p95", "n_perm", "emp_p"
    )
}
