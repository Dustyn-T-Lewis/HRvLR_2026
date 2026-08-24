# What |rho| does a protein reach against the blood index by chance?
#
# blood_cor_max gates removal, so the threshold needs a null rather than a round
# number. Permuting the index across samples breaks the protein-to-blood
# association while leaving both marginal distributions and the matrix's own
# missingness alone, so the pooled |rho| across permutations is the distribution
# a single protein draws from when it tracks nothing.
#
# filter_config.R quotes this quantile. Recompute it with
# 01_Filtering/a_script/00_blood_cor_null.R whenever the input matrix changes -
# a cut inherited from a matrix that no longer exists is one nobody can defend.

pacman::p_load(stats)

blood_cor_null_quantile <- function(log_int, blood_index, n_perm = 300,
                                    prob = 0.999, seed = 42) {
  set.seed(seed)
  null_rho <- vapply(seq_len(n_perm), function(i) {
    permuted <- blood_index[sample.int(length(blood_index))]
    abs(as.vector(suppressWarnings(stats::cor(
      t(log_int), permuted,
      method = "spearman", use = "pairwise.complete.obs"
    ))))
  }, numeric(nrow(log_int)))
  stats::quantile(null_rho, prob, na.rm = TRUE, names = FALSE)
}
