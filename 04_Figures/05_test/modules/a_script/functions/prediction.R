pacman::p_load(pROC, withr)

module_auc <- function(y, x) {
  ok <- is.finite(x) & !is.na(y)
  a <- as.numeric(pROC::auc(pROC::roc(y[ok], x[ok], quiet = TRUE, direction = "auto")))
  max(a, 1 - a)
}

perm_p_unpaired <- function(y, x, n_perm = 1000, seed = 42) {
  obs <- module_auc(y, x)
  withr::with_seed(seed, {
    null <- vapply(seq_len(n_perm), function(i) module_auc(sample(y), x), numeric(1))
  })
  (sum(null >= obs) + 1) / (n_perm + 1)
}
