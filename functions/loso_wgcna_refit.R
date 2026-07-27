# Refit the network inside each leave-one-subject-out fold to measure how much
# module definition borrows from the subject it is asked to predict. The
# cohort eigengenes are computed once over all 16 subjects, so a module used
# to predict a held-out subject was partly defined by that subject. Rebuilding
# per fold removes that, and the drop between the two is the size of the
# circularity.
#
# The held-out eigengene comes from a PC1 fitted on training columns only;
# WGCNA::moduleEigengenes over the full matrix would see the held-out rows.

pacman::p_load(here, dplyr, tibble)
source(here("functions", "shared_wgcna.R"))

MIN_JACCARD <- 0.5

fit_pc1 <- function(x) {
  if (ncol(x) < 2) {
    stop("fit_pc1 needs at least two proteins")
  }
  center <- colMeans(x)
  scale <- apply(x, 2, stats::sd)
  scale[scale == 0] <- 1
  z <- scale(x, center = center, scale = scale)
  v <- svd(z)$v[, 1]
  if (stats::cor(as.numeric(z %*% v), rowMeans(z)) < 0) {
    v <- -v
  }
  list(center = center, scale = scale, rotation = v)
}

project_pc1 <- function(fit, x) {
  z <- scale(x, center = fit$center, scale = fit$scale)
  as.numeric(z %*% fit$rotation)
}

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

match_modules <- function(mod_train, mod_full) {
  full_lvls <- setdiff(unique(mod_full), "grey")
  train_lvls <- setdiff(unique(mod_train), "grey")

  rows <- lapply(full_lvls, function(f) {
    a <- names(mod_full)[mod_full == f]
    j <- vapply(train_lvls, function(t) {
      jaccard(a, names(mod_train)[mod_train == t])
    }, numeric(1))
    data.frame(
      full = f, train = train_lvls[[which.max(j)]], jaccard = max(j),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

# One fold: rebuild the network without `hold`, match its modules back to the
# cohort's, and project the held-out samples onto the rebuilt module.
refit_fold <- function(abund, subject, hold, modules_full,
                       min_jaccard = MIN_JACCARD) {
  train <- subject != hold
  net <- fit_modules(abund[, train, drop = FALSE], subject[train])
  mm <- match_modules(net$colors, modules_full)

  me_hold <- lapply(mm$full, function(f) {
    row <- mm[mm$full == f, ]
    if (row$jaccard < min_jaccard) {
      return(NULL)
    }
    prot <- names(net$colors)[net$colors == row$train]
    if (length(prot) < 2) {
      return(NULL)
    }
    fit <- fit_pc1(t(abund[prot, train, drop = FALSE]))
    project_pc1(fit, t(abund[prot, !train, drop = FALSE]))
  })
  names(me_hold) <- mm$full

  list(me_hold = Filter(Negate(is.null), me_hold), jaccard = mm)
}
