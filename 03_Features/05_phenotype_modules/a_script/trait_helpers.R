# Label-free feature-trait association machinery shared by the three
# 05_phenotype_modules scripts. The design is fixed before results are
# read: Spearman rho per feature x phenotype x timepoint, a
# subject-permutation null with one shared permutation index so replicate
# b is the same shuffled cohort in every timepoint, and a consistency
# call that requires the same sign in all timepoints and empirical
# p < 0.05 in at least two. Timepoints share subjects, so consistency is
# stability of an association, not independent replication.

perm_index <- function(n_subj, n_perm = 1000, seed = 42) {
  set.seed(seed)
  replicate(n_perm, sample.int(n_subj))
}

# feature_mat: features x subjects, columns aligned to pheno. Returns
# observed rho, the empirical two-sided p, and the null rho matrix the
# consistency null reuses. Pairwise-complete handles the one d_1rm_ext
# NA; the permutation moves the NA with the value, preserving its count.
cor_scan <- function(feature_mat, pheno, perm_idx) {
  spear <- function(y) {
    suppressWarnings(cor(t(feature_mat), y,
      method = "spearman", use = "pairwise.complete.obs"
    ))
  }
  obs <- as.vector(spear(pheno))
  null_rho <- spear(matrix(pheno[perm_idx], nrow = length(pheno)))
  list(
    obs = tibble::tibble(
      feature = rownames(feature_mat),
      rho = obs,
      emp_p = (rowSums(abs(null_rho) >= abs(obs)) + 1) / (ncol(null_rho) + 1)
    ),
    null_rho = null_rho
  )
}

# Per-replicate empirical p of each null rho against its own row, so the
# consistency criterion can be applied to every permuted cohort exactly
# as it is to the observed one.
null_emp_p <- function(null_rho) {
  t(apply(-abs(null_rho), 1, rank, ties.method = "max")) / ncol(null_rho)
}

consistent_call <- function(rho_mat, p_mat, alpha = 0.05, min_sig = 2L) {
  same_sign <- apply(sign(rho_mat), 1, function(s) {
    all(s != 0) && length(unique(s)) == 1
  })
  same_sign & rowSums(p_mat < alpha) >= min_sig
}

# Observed consistency count against the count in each permuted cohort.
# scans: list over timepoints from cor_scan for one phenotype.
consistency_scan <- function(scans, alpha = 0.05, min_sig = 2L) {
  rho_obs <- vapply(scans, function(s) s$obs$rho, scans[[1]]$obs$rho)
  p_obs <- vapply(scans, function(s) s$obs$emp_p, scans[[1]]$obs$emp_p)
  obs_call <- consistent_call(rho_obs, p_obs, alpha, min_sig)

  null_p <- lapply(scans, function(s) null_emp_p(s$null_rho))
  n_perm <- ncol(scans[[1]]$null_rho)
  null_counts <- vapply(seq_len(n_perm), function(b) {
    rho_b <- vapply(scans, function(s) s$null_rho[, b], rho_obs[, 1])
    p_b <- vapply(null_p, function(p) p[, b], p_obs[, 1])
    sum(consistent_call(rho_b, p_b, alpha, min_sig))
  }, integer(1))

  list(
    consistent = scans[[1]]$obs$feature[obs_call],
    n_observed = sum(obs_call),
    null_median = stats::median(null_counts),
    emp_p = (sum(null_counts >= sum(obs_call)) + 1) / (n_perm + 1)
  )
}

# Silhouette-maximising k for hierarchical clustering on a distance.
choose_k <- function(d, hc, ks = 2:10) {
  sil <- vapply(ks, function(k) {
    mean(cluster::silhouette(stats::cutree(hc, k), d)[, "sil_width"])
  }, numeric(1))
  ks[which.max(sil)]
}
