pacman::p_load(dplyr, tibble, purrr)

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

# Scaled PC1 loading, sign-oriented so the score tracks the row mean.
fit_pc1 <- function(train_mat) {
  ctr <- colMeans(train_mat)
  scl <- apply(train_mat, 2, sd)
  scl[scl == 0 | !is.finite(scl)] <- 1
  z <- sweep(sweep(train_mat, 2, ctr, "-"), 2, scl, "/")
  v <- svd(z, nu = 0, nv = 1)$v[, 1]
  score <- as.numeric(z %*% v)
  if (cor(score, rowMeans(z)) < 0) v <- -v
  list(center = ctr, scale = scl, loading = v)
}

project_pc1 <- function(fit, new_mat) {
  z <- sweep(sweep(new_mat, 2, fit$center, "-"), 2, fit$scale, "/")
  as.numeric(z %*% fit$loading)
}

match_modules <- function(colors_full, colors_train) {
  full_levels <- setdiff(unique(colors_full), "grey")
  train_levels <- setdiff(unique(colors_train), "grey")
  purrr::map_dfr(full_levels, function(f) {
    prot_f <- names(colors_full)[colors_full == f]
    j <- vapply(train_levels, function(t) {
      jaccard(prot_f, names(colors_train)[colors_train == t])
    }, numeric(1))
    tibble(full = f, train = train_levels[which.max(j)], jaccard = max(j))
  })
}
