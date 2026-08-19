# Shared assembly and statistics for the mCSA axis stage (PREREG.md). Every
# script here reads the proteome through pilot_data(), so the complete-case
# gate and the blood-index join stay the ones the galamm pilot used.

pacman::p_load(here, dplyr, tibble, readr, limma, withr)

CSA_ITEMS <- c("d_fcsa_I", "d_fcsa_II", "d_mcsa")
MCSA_CONFIGS <- c("T1", "T2", "T3", "delta")
DISCORDANT_SUBJECT <- "LR_S14"

phenotype_table <- function() {
  read_csv(
    here("00_input", "c_data", "phenotype.csv"),
    show_col_types = FALSE
  )
}

# The axis the two fibre measures share (r = 0.91), on a common scale.
fibre_axis <- function(pheno) {
  as.numeric(rowMeans(scale(pheno[, c("d_fcsa_I", "d_fcsa_II")])))
}

# What each CSA measure contributes to the composite the arms were split on,
# and what the fibre pair alone recovers without the whole-muscle term.
composite_weights <- function(pheno) {
  z <- as.data.frame(scale(pheno[, CSA_ITEMS]))
  z$comp <- pheno$comp_hypertrophy
  full <- stats::lm(comp ~ d_fcsa_I + d_fcsa_II + d_mcsa, z)
  fibre <- stats::lm(comp ~ d_fcsa_I + d_fcsa_II, z)
  tibble(
    item = CSA_ITEMS,
    weight = unname(stats::coef(full)[CSA_ITEMS]),
    r2_full = summary(full)$r.squared,
    r2_without_mcsa = summary(fibre)$r.squared
  )
}

# Probability a random HR outranks a random LR: the rank statistic behind the
# Wilcoxon test, on the scale a reader can compare to 0.5. d_1rm_ext has one NA,
# so the subject is dropped rather than propagating through the comparison.
arm_auc <- function(value, arm) {
  ok <- !is.na(value)
  hr <- value[ok & arm == "HR"]
  lr <- value[ok & arm == "LR"]
  mean(outer(hr, lr, ">") + 0.5 * outer(hr, lr, "=="))
}

# One design per config. T1/T2/T3 take that timepoint's samples; delta takes
# the per-subject T2 - T1 difference in both abundance and blood index, so the
# covariate is differenced on the same footing as the response.
config_design <- function(mat, meta, pheno, config, drop = NULL) {
  keep <- !meta$subject %in% drop
  meta <- meta[keep, , drop = FALSE]
  mat <- mat[, keep, drop = FALSE]

  if (config == "delta") {
    subject <- intersect(
      meta$subject[meta$timepoint == "T1"], meta$subject[meta$timepoint == "T2"]
    )
    at <- function(tp) {
      match(paste(subject, tp), paste(meta$subject, meta$timepoint))
    }
    late <- at("T2")
    early <- at("T1")
    y <- mat[, late, drop = FALSE] - mat[, early, drop = FALSE]
    blood <- meta$blood_index[late] - meta$blood_index[early]
  } else {
    i <- which(meta$timepoint == config)
    subject <- meta$subject[i]
    y <- mat[, i, drop = FALSE]
    blood <- meta$blood_index[i]
  }

  list(
    y = y, blood = blood, subject = subject,
    mcsa = pheno$d_mcsa[match(subject, pheno$subject)]
  )
}

# limma with a continuous predictor: one moderated t per protein for the d_mcsa
# slope with the blood index partialled out, BH across the 931 within a config.
# One sample per subject here, so there are no repeated measures to block on.
mcsa_scan <- function(design) {
  x <- stats::model.matrix(
    ~ mcsa + blood, data.frame(mcsa = design$mcsa, blood = design$blood)
  )
  fit <- eBayes(lmFit(design$y, x))
  topTable(fit, coef = "mcsa", number = Inf, sort.by = "none") |>
    as_tibble(rownames = "feature") |>
    transmute(
      feature = .data$feature, slope = .data$logFC, t = .data$t,
      p = .data$P.Value, bh = .data$adj.P.Val
    )
}

# What a survivor rests on: how far the slope's t moves when any one subject
# leaves, whether it holds inside each arm, and whether the arm label explains
# it. d_mcsa separates the arms at AUC 0.83, so a continuous hit could be the
# group contrast in another coat.
survivor_checks <- function(design, feature) {
  y <- design$y[feature, ]
  slope_t <- function(i) {
    stats::coef(summary(stats::lm(
      y[-i] ~ design$mcsa[-i] + design$blood[-i]
    )))[2, 3]
  }
  loo <- vapply(seq_along(y), slope_t, numeric(1))
  arm <- as.numeric(startsWith(design$subject, "HR"))
  both <- stats::coef(summary(
    stats::lm(y ~ design$mcsa + design$blood + arm)
  ))
  arm_rho <- function(a) {
    stats::cor(y[arm == a], design$mcsa[arm == a], method = "spearman")
  }
  tibble(
    feature = feature,
    loo_t_min = min(loo), loo_t_max = max(loo),
    rho_hr = arm_rho(1), rho_lr = arm_rho(0),
    slope_with_arm = both[2, 1], p_with_arm = both[2, 4],
    p_arm_term = both[4, 4]
  )
}

# Shuffle d_mcsa across subjects and refit. Subject is the unit because the
# phenotype is a subject property; abundances and blood index stay put.
mcsa_permutation <- function(design, observed, n_perm = 200, seed = 42) {
  null_t <- with_seed(seed, vapply(seq_len(n_perm), function(i) {
    permuted <- design
    permuted$mcsa <- sample(design$mcsa)
    scan <- mcsa_scan(permuted)
    abs(scan$t[match(observed$feature, scan$feature)])
  }, numeric(nrow(observed))))

  # vapply drops to a vector when only one feature survived; rowSums needs both.
  n_ge <- rowSums(matrix(null_t, nrow = nrow(observed)) >= abs(observed$t))
  tibble(
    feature = observed$feature, obs_t = observed$t, n_ge = n_ge,
    n_perm = n_perm, emp_p = (n_ge + 1) / (n_perm + 1)
  )
}
