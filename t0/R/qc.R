# Quality control on the raw modeling inputs. Each returns a tibble the report renders; outlier
# checks write flags rather than dropping samples, so the decision stays visible.

qc_missingness <- function(proteins, meta) {
  expr <- proteins$expr
  per_protein <- tibble::tibble(
    protein = rownames(expr),
    detected = matrixStats::rowSums2(!is.na(expr)),
    frac = detected / ncol(expr)
  )
  per_sample <- tibble::tibble(
    sample = colnames(expr),
    detected = matrixStats::colSums2(!is.na(expr)),
    frac = detected / nrow(expr)
  ) |>
    dplyr::left_join(dplyr::select(meta, sample, subject, group, timepoint), by = "sample")
  list(per_protein = per_protein, per_sample = per_sample)
}

qc_distributions <- function(proteins, meta) {
  expr <- proteins$expr
  tibble::tibble(
    sample = colnames(expr),
    median = matrixStats::colMedians(expr, na.rm = TRUE),
    iqr = matrixStats::colIQRs(expr, na.rm = TRUE),
    cv = matrixStats::colSds(expr, na.rm = TRUE) / matrixStats::colMeans2(expr, na.rm = TRUE)
  ) |>
    dplyr::left_join(dplyr::select(meta, sample, group, timepoint), by = "sample")
}

qc_correlation <- function(proteins) {
  expr <- proteins$expr
  complete <- expr[matrixStats::rowSums2(is.na(expr)) == 0, , drop = FALSE]
  stats::cor(complete, method = "pearson")
}

qc_pairing <- function(meta, cfg) {
  meta |>
    dplyr::count(subject, timepoint) |>
    tidyr::pivot_wider(names_from = timepoint, values_from = n, values_fill = 0) |>
    dplyr::mutate(
      complete_triad = rowSums(dplyr::across(dplyr::all_of(cfg$design$timepoints)) > 0) ==
        length(cfg$design$timepoints)
    )
}

# Two independent outlier signals on the complete-case matrix: how far a sample's median correlation
# to the others sits below the group, and its Mahalanobis distance in the leading PCs.
qc_outliers <- function(proteins, cfg) {
  cor_mat <- qc_correlation(proteins)
  med_cor <- apply(cor_mat, 2, function(x) stats::median(x[x < 1]))
  cor_z <- (med_cor - stats::median(med_cor)) / stats::mad(med_cor)

  complete <- proteins$expr[matrixStats::rowSums2(is.na(proteins$expr)) == 0, , drop = FALSE]
  pcs <- stats::prcomp(t(complete), scale. = TRUE)$x[, 1:3]
  mahal <- stats::mahalanobis(pcs, colMeans(pcs), stats::cov(pcs))

  tibble::tibble(
    sample = colnames(cor_mat),
    median_cor = med_cor,
    cor_z = cor_z,
    mahal = mahal,
    flag_cor = cor_z < cfg$thresholds$outlier_cor_z,
    flag_pca = mahal > stats::qchisq(1 - cfg$thresholds$outlier_mahal_p, df = 3),
    flagged = flag_cor | flag_pca
  )
}
