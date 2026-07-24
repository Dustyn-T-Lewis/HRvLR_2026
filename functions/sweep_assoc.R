# In-sample association for the screening sweep. Each cell is one feature space,
# one timepoint config, one method, one outcome: a per-feature test of whether
# the feature moves with the outcome across the 16 subjects, with BH within the
# cell. This is descriptive; the prediction roots adjudicate generalization.
# Feature matrices are the same subject-level bundle the prediction roots use,
# so proteins carry the same cohort-imputed optimism label throughout.

pacman::p_load(dplyr, tibble, limma)

# Ordinary simple regression of the outcome on each feature, vectorised across
# features. Returns the slope, its t, the two-sided p, and BH within the cell.
assoc_lm <- function(x_mat, y) {
  n <- length(y)
  yc <- y - mean(y)
  xc <- sweep(x_mat, 2, colMeans(x_mat))
  sxx <- colSums(xc^2)
  sxy <- colSums(xc * yc)
  slope <- sxy / sxx
  ss_res <- sum(yc^2) - slope * sxy
  se <- sqrt((ss_res / (n - 2)) / sxx)
  t <- unname(slope / se)
  tibble(
    feature = colnames(x_mat), effect = unname(slope), t = t,
    p = 2 * stats::pt(-abs(t), df = n - 2)
  ) |>
    mutate(bh = stats::p.adjust(.data$p, "BH"))
}

# Spearman rank association: rho and the t-approximation p per feature.
assoc_spearman <- function(x_mat, y) {
  n <- length(y)
  rho <- as.numeric(stats::cor(apply(x_mat, 2, rank), rank(y)))
  denom <- pmax(1 - rho^2, .Machine$double.eps)
  t <- rho * sqrt((n - 2) / denom)
  tibble(
    feature = colnames(x_mat), effect = rho, t = t,
    p = 2 * stats::pt(-abs(t), df = n - 2)
  ) |>
    mutate(bh = stats::p.adjust(.data$p, "BH"))
}

# Moderated association (limma eBayes) of each feature on a design. coef 2 is
# the effect of interest (the slope for a continuous outcome, the HR-LR
# difference for a group design).
assoc_limma_fit <- function(x_mat, design) {
  fit <- limma::eBayes(limma::lmFit(t(x_mat), design))
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tibble(
    feature = rownames(tt), effect = tt$logFC, t = tt$t,
    p = tt$P.Value, bh = tt$adj.P.Val
  )
}

assoc_limma_cont <- function(x_mat, y) {
  assoc_limma_fit(x_mat, stats::model.matrix(~y))
}

assoc_limma_group <- function(x_mat, y) {
  g <- factor(ifelse(y == 1, "HR", "LR"), levels = c("LR", "HR"))
  assoc_limma_fit(x_mat, stats::model.matrix(~g))
}

# Wilcoxon rank-sum HR-vs-LR difference per feature, effect as the median shift.
assoc_wilcox_group <- function(x_mat, y) {
  hr <- y == 1
  stat <- vapply(seq_len(ncol(x_mat)), function(j) {
    col <- x_mat[, j]
    w <- suppressWarnings(stats::wilcox.test(col[hr], col[!hr]))
    c(effect = stats::median(col[hr]) - stats::median(col[!hr]), p = w$p.value)
  }, numeric(2))
  tibble(
    feature = colnames(x_mat), effect = stat["effect", ], t = NA_real_,
    p = stat["p", ]
  ) |>
    mutate(bh = stats::p.adjust(.data$p, "BH"))
}

# One association cell. is_group switches between the continuous methods
# (limma / lm / spearman) and the group HR-vs-LR methods (limma, wilcoxon).
run_assoc_cell <- function(x_mat, y, method, is_group) {
  if (is_group) {
    switch(method,
      limma = assoc_limma_group(x_mat, y),
      wilcoxon = assoc_wilcox_group(x_mat, y),
      stop("unknown group method: ", method)
    )
  } else {
    switch(method,
      limma = assoc_limma_cont(x_mat, y),
      lm = assoc_lm(x_mat, y),
      spearman = assoc_spearman(x_mat, y),
      stop("unknown continuous method: ", method)
    )
  }
}
