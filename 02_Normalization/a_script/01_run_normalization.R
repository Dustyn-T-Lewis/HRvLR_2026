#!/usr/bin/env Rscript
# HRvLR Stage 02: cycloess normalization of the filtered, non-imputed DAList.
# limma handles per-protein NAs downstream, so the canonical Stage 03 input stays
# non-imputed; the imputation/ arms add exploratory imputed DALists alongside.
# cycloess is limma::normalizeCyclicLoess(method = "fast"), 3 iterations. The span is not
# 0.7: limma defaults adaptive.span = TRUE, which overrides the 0.7 formal with
# chooseLowessSpan(nrow) = 0.5075 at 1920 proteins. It moves if the protein count moves.

pacman::p_load(proteoDA, here, readr, dplyr, stringr, tibble)
source(here("04_Figures", "functions", "shared_pca.R"))
source(here("04_Figures", "functions", "shared_utils.R"))

data_dir <- here("02_Normalization", "c_data")
report_dir <- here("02_Normalization", "b_reports")
clear_dir(data_dir)
clear_dir(report_dir)

# Load filtered DAList

dal <- readRDS(here("01_Filtering", "c_data", "DAList_filtered.rds"))
cat(sprintf("Loaded filtered DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# Normalize (cycloess) + QC

write_norm_report(dal,
  grouping_column = "Group_Time", output_dir = report_dir,
  filename = "norm_comparison.pdf", overwrite = TRUE
)
write_qc_report(dal,
  color_column = "Group_Time", output_dir = report_dir,
  filename = "qc_pre.pdf", overwrite = TRUE
)
dal <- normalize_data(dal, norm_method = "cycloess")
write_qc_report(dal,
  color_column = "Group_Time", output_dir = report_dir,
  filename = "qc_post.pdf", overwrite = TRUE
)
cat(sprintf("Normalized (cycloess): %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# Export

export_df <- bind_cols(
  as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
  as_tibble(dal$data)
)
write_csv(export_df, file.path(data_dir, "normalized.csv"))
saveRDS(dal, file.path(data_dir, "DAList_normalized.rds"))

# Combined report intermediates (filtering + normalization)

subj_var <- dal$metadata |>
  mutate(
    iqr = apply(dal$data[, Col_ID], 2, IQR, na.rm = TRUE),
    Subject_ID = str_remove(Col_ID, "_T[123]$")
  ) |>
  select(Col_ID, Subject_ID, Group, Timepoint, Group_Time, iqr)

log_dat <- dal$data
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
  dal_nrow = nrow(dal$data), dal_ncol = ncol(dal$data)
))
saveRDS(intermediates, file.path(data_dir, "00_report_intermediates.rds"))

if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Done: %d proteins x %d samples -> %s/\n", nrow(dal$data), ncol(dal$data), data_dir))
