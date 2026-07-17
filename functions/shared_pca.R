# Shared PCA helper for the filtering and normalization stages: median-fills gaps,
# optionally log2-transforms, and returns the prcomp fit, the PC1-3 scores joined
# to metadata, and the top-3 variance-explained percentages.

run_pca <- function(mat, metadata, log_transform = TRUE) {
  for (j in seq_len(ncol(mat))) {
    mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  }
  if (log_transform) mat <- log2(mat + 1)
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  ve <- round(summary(pca)$importance[2, 1:3] * 100, 1)
  pc <- as.data.frame(pca$x[, 1:3]) |>
    mutate(Col_ID = rownames(pca$x)) |>
    left_join(metadata, by = "Col_ID")
  list(pca = pca, scores = pc, var_exp = ve)
}
