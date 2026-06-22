#!/usr/bin/env Rscript
# HRvLR Stage 02: cycloess normalization of the filtered, non-imputed DAList.
# limma handles per-protein NAs downstream, so the canonical Stage 03 input stays
# non-imputed; the imputation/ arms add exploratory imputed DALists alongside.
# cycloess uses limma's defaults (span = 0.7, adaptive.span = FALSE), matching YvO.

pacman::p_load(proteoDA, here, readr, dplyr, tidyr, stringr, tibble)
set.seed(42)

data_dir <- here("02_Normalization", "c_data")
report_dir <- here("02_Normalization", "b_reports")
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}
clear_dir(data_dir)
clear_dir(report_dir)

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

#### Load filtered DAList ####

dal <- readRDS(here("01_Filtering", "c_data", "DAList_filtered.rds"))
cat(sprintf("Loaded filtered DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

#### Normalize (cycloess) + QC ####

write_norm_report(dal,
  grouping_column = "Group_Time", output_dir = report_dir,
  filename = "norm_comparison.pdf", overwrite = TRUE
)
write_qc_report(dal,
  color_column = "Group_Time", output_dir = report_dir,
  filename = "qc_pre.pdf", overwrite = TRUE
)
dal_pre <- dal
dal <- normalize_data(dal, norm_method = "cycloess")
write_qc_report(dal,
  color_column = "Group_Time", output_dir = report_dir,
  filename = "qc_post.pdf", overwrite = TRUE
)
cat(sprintf("Normalized (cycloess): %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

#### Normalization method ranking (PRONE-style) ####

norm_metric <- function(mat, groups, metric) {
  grp_list <- split(seq_len(ncol(mat)), groups)
  if (metric == "cor") {
    vals <- unlist(lapply(grp_list, function(idx) {
      sub <- mat[, idx, drop = FALSE]
      if (ncol(sub) < 2) {
        return(numeric(0))
      }
      cm <- cor(sub, use = "pairwise.complete.obs")
      cm[lower.tri(cm)]
    }))
    return(mean(vals, na.rm = TRUE))
  }
  vals <- unlist(lapply(grp_list, function(idx) {
    sub <- mat[, idx, drop = FALSE]
    apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) {
        return(NA_real_)
      }
      if (metric == "cv") sd(x) / abs(mean(x)) else mad(x, constant = 1)
    })
  }))
  median(vals, na.rm = TRUE)
}

dal_pre$metadata$group <- factor(dal_pre$metadata$Group_Time)
methods <- c("log2", "median", "mean", "vsn", "quantile", "cycloess", "rlr", "gi")
norm_scores <- lapply(methods, function(m) {
  dal_n <- tryCatch(normalize_data(dal_pre, norm_method = m), error = function(e) NULL)
  if (is.null(dal_n)) {
    return(NULL)
  }
  mat <- as.matrix(dal_n$data)
  grps <- dal_n$metadata$group
  tibble(
    method = m,
    PCV = round(norm_metric(mat, grps, "cv"), 4),
    PMAD = round(norm_metric(mat, grps, "mad"), 4),
    COR = round(norm_metric(mat, grps, "cor"), 4)
  )
}) |>
  bind_rows() |>
  mutate(
    PCV_rank = rank(PCV), PMAD_rank = rank(PMAD), COR_rank = rank(-COR),
    composite = round((PCV_rank + PMAD_rank + COR_rank) / 3, 2)
  ) |>
  arrange(composite)
write_csv(norm_scores, file.path(data_dir, "norm_quality_scores.csv"))
cat("Norm ranking:\n")
print(norm_scores |> select(method, PCV, PMAD, COR, composite))

#### Export ####

export_df <- bind_cols(
  as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
  as_tibble(dal$data)
)
write_csv(export_df, file.path(data_dir, "normalized.csv"))
saveRDS(dal, file.path(data_dir, "DAList_normalized.rds"))

#### Combined report intermediates (filtering + normalization) ####

subj_var <- dal$metadata |>
  mutate(
    iqr = apply(log2(dal$data[, Col_ID] + 1), 2, IQR, na.rm = TRUE),
    Subject_ID = str_remove(Col_ID, "_T[123]$")
  ) |>
  select(Col_ID, Subject_ID, Group, Timepoint, Group_Time, iqr)

log_dat <- log2(dal$data + 1)
grp_vec <- dal$metadata$Group_Time[match(colnames(log_dat), dal$metadata$Col_ID)]
eta2_vals <- apply(log_dat, 1, function(x) {
  ok <- !is.na(x)
  if (sum(ok) < 4) {
    return(NA_real_)
  }
  xk <- x[ok]
  gk <- grp_vec[ok]
  ss_b <- sum(tapply(xk, gk, length) * (tapply(xk, gk, mean) - mean(xk))^2)
  ss_t <- sum((xk - mean(xk))^2)
  if (ss_t > 0) ss_b / ss_t else NA_real_
})
pca_post <- run_pca(dal$data, dal$metadata, log_transform = FALSE)

filt <- readRDS(here("01_Filtering", "c_data", "filtering_intermediates.rds"))
intermediates <- c(filt, list(
  pca_post = pca_post, subj_var = subj_var, eta2_vals = eta2_vals,
  norm_scores = norm_scores, dal_nrow = nrow(dal$data), dal_ncol = ncol(dal$data)
))
saveRDS(intermediates, file.path(data_dir, "00_report_intermediates.rds"))

if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Done: %d proteins x %d samples -> %s/\n", nrow(dal$data), ncol(dal$data), data_dir))
