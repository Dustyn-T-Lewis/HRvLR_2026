#!/usr/bin/env Rscript
# Canonical HRvLR Stage 02 imputation
# Direct YvO transfer: 3-method MAR/MNAR consensus -> missForest
#
# Primary DEP remains non-imputed. This stage exists for QC, sensitivity, and
# complete-matrix downstream tasks.
#
# Outputs:
#   02_Imputation/c_data/00_report_intermediates.rds
#   02_Imputation/c_data/01_imputed.csv
#   02_Imputation/c_data/01_DAList_imputed.rds
#   02_Imputation/c_data/02_mar_mnar_classification.csv
#   02_Imputation/c_data/07_imputation_mask.csv
#   02_Imputation/c_data/08_mnar_imputation_audit.csv
#   02_Imputation/c_data/09_imputation_summary.txt
#   02_Imputation/c_data/sessionInfo.txt

library(missForest)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)

set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

cfg <- list(
  norm_csv = "01_normalization/c_data/02_normalized.csv",
  norm_rds = "01_normalization/c_data/03_DAList_normalized.rds",
  data_dir = "02_Imputation/c_data",
  miss_unreliable = 50
)

dir.create(cfg$data_dir, recursive = TRUE, showWarnings = FALSE)

PAL_GT <- c(
  HR_T1 = scales::alpha("#D6604D", 0.4),
  HR_T2 = "#D6604D",
  HR_T3 = scales::alpha("#B2182B", 0.8),
  LR_T1 = scales::alpha("#4393C3", 0.4),
  LR_T2 = "#4393C3",
  LR_T3 = scales::alpha("#2166AC", 0.8)
)
PAL_MAR <- c(MAR = "#4393C3", MNAR = "#D6604D")
PAL_CLASS <- c(Complete = "#4DAF4A", MAR = "#4393C3", MNAR = "#D6604D")

required_meta_cols <- c("Col_ID", "Subject_ID", "Group", "Timepoint", "Group_Time")

df <- read_csv(cfg$norm_csv, show_col_types = FALSE)
ann <- df |>
  select(uniprot_id, gene, protein, description)

feature_id <- if (anyDuplicated(df$gene)) df$uniprot_id else df$gene
ann <- ann |>
  mutate(feature_id = feature_id)

mat <- as.matrix(df[, -(1:4)])
rownames(mat) <- feature_id

dal_norm <- readRDS(cfg$norm_rds)
missing_meta_cols <- setdiff(required_meta_cols, names(dal_norm$metadata))
if (length(missing_meta_cols) > 0) {
  stop(sprintf(
    "Missing required metadata columns in normalized DAList: %s",
    paste(missing_meta_cols, collapse = ", ")
  ))
}

meta <- as_tibble(dal_norm$metadata) |>
  select(all_of(required_meta_cols))

stopifnot("Sample mismatch" = setequal(meta$Col_ID, colnames(mat)))
mat <- mat[, meta$Col_ID, drop = FALSE]

meta <- meta |>
  mutate(
    Group = factor(Group, levels = c("HR", "LR")),
    Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3")),
    Group_Time = factor(Group_Time,
      levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3"))
  )

if (any(is.na(meta$Subject_ID))) {
  stop("Subject_ID contains missing values; repeated-measures imputation audit cannot proceed.")
}

cat(sprintf("Loaded: %d proteins x %d samples\n", nrow(mat), ncol(mat)))

prot_miss <- rowSums(is.na(mat))
prot_pct <- prot_miss / ncol(mat) * 100
obs_means <- rowMeans(mat, na.rm = TRUE)
pct_miss <- round(sum(is.na(mat)) / length(mat) * 100, 2)

miss_by_group <- sapply(unique(as.character(meta$Group_Time)), function(gt) {
  cols <- meta$Col_ID[meta$Group_Time == gt]
  rowSums(is.na(mat[, cols, drop = FALSE])) / length(cols) * 100
})

has_na <- which(prot_miss > 0 & prot_miss < ncol(mat))
inc_mean <- obs_means[has_na]
inc_pct <- prot_pct[has_na]

# Classifier 1: k-means on observed intensity and missingness
set.seed(42)
km <- kmeans(scale(cbind(inc_mean, inc_pct)), centers = 2, nstart = 25)
km_mnar <- km$cluster == which.min(tapply(inc_mean, km$cluster, mean))

# Classifier 2: global logistic P(missing | feature mean intensity)
lr_df <- data.frame(
  is_miss = as.integer(is.na(as.vector(mat))),
  intensity = rep(obs_means, ncol(mat))
)
lr_fit <- glm(is_miss ~ intensity, data = lr_df, family = binomial)
lr_pred <- predict(lr_fit, newdata = data.frame(intensity = inc_mean), type = "response")
lr_mnar <- lr_pred > median(lr_pred)

# Classifier 3: left-tail proximity
global_q25 <- quantile(mat, 0.25, na.rm = TRUE)
tail_frac <- vapply(
  has_na,
  function(i) mean(mat[i, !is.na(mat[i, ])] < global_q25),
  numeric(1)
)
lt_mnar <- (tail_frac * inc_pct / 100) > median(tail_frac * inc_pct / 100)

votes <- as.integer(km_mnar) + as.integer(lr_mnar) + as.integer(lt_mnar)

miss_class <- tibble(
  uniprot_id = ann$uniprot_id,
  gene = ann$gene,
  protein = ann$protein,
  description = ann$description,
  feature_id = ann$feature_id,
  n_miss = prot_miss,
  pct_miss = prot_pct,
  mean_intensity = obs_means,
  vote_kmeans = NA_integer_,
  vote_logistic = NA_integer_,
  vote_lefttail = NA_integer_,
  n_mnar_votes = NA_integer_
)

miss_class$vote_kmeans[has_na] <- as.integer(km_mnar)
miss_class$vote_logistic[has_na] <- as.integer(lr_mnar)
miss_class$vote_lefttail[has_na] <- as.integer(lt_mnar)
miss_class$n_mnar_votes[has_na] <- votes

miss_class <- miss_class |>
  mutate(
    classification = case_when(
      n_miss == 0 ~ "Complete",
      n_miss >= ncol(mat) ~ "MNAR",
      n_mnar_votes >= 2 ~ "MNAR",
      TRUE ~ "MAR"
    ),
    imputation_reliable = classification == "Complete" | pct_miss < cfg$miss_unreliable
  )

mnar_ids <- miss_class$feature_id[miss_class$classification == "MNAR"]

group_miss_pval <- setNames(rep(NA_real_, nrow(miss_class)), miss_class$feature_id)
for (fid in mnar_ids) {
  ct <- sapply(unique(as.character(meta$Group_Time)), function(gt) {
    cols <- meta$Col_ID[meta$Group_Time == gt]
    c(missing = sum(is.na(mat[fid, cols])), observed = sum(!is.na(mat[fid, cols])))
  })
  group_miss_pval[fid] <- tryCatch(
    fisher.test(ct, simulate.p.value = TRUE, B = 2000)$p.value,
    error = function(e) NA_real_
  )
}
miss_class$group_miss_pval <- group_miss_pval[miss_class$feature_id]

n_mar <- sum(miss_class$classification == "MAR")
n_mnar <- sum(miss_class$classification == "MNAR")
n_complete <- sum(miss_class$classification == "Complete")
mar_vals <- sum(miss_class$n_miss[miss_class$classification == "MAR"])
mnar_vals <- sum(miss_class$n_miss[miss_class$classification == "MNAR"])
total_vals <- mar_vals + mnar_vals
pct_mar_vals <- if (total_vals > 0) round(mar_vals / total_vals * 100, 1) else NA_real_

cat(sprintf(
  "Classification: MAR %d | MNAR %d | Complete %d\n",
  n_mar, n_mnar, n_complete
))

cat("Imputing with missForest...\n")
gene_order <- order(rownames(mat))
mat <- mat[gene_order, , drop = FALSE]
ann <- ann[gene_order, , drop = FALSE]
miss_class <- miss_class[gene_order, , drop = FALSE]

set.seed(42)
mf <- missForest::missForest(t(mat), maxiter = 10, ntree = 100, verbose = TRUE)
mat_imp <- t(mf$ximp)
rownames(mat_imp) <- rownames(mat)
colnames(mat_imp) <- colnames(mat)
stopifnot("Imputed matrix still contains NA values" = sum(is.na(mat_imp)) == 0)

oob <- as.numeric(mf$OOBerror[1])
cat(sprintf("OOB error: %.4f\n", oob))

was_na <- is.na(mat)
mnar_idx <- miss_class$feature_id[miss_class$classification == "MNAR"]

mnar_audit <- tibble(
  uniprot_id = miss_class$uniprot_id[match(mnar_idx, miss_class$feature_id)],
  gene = miss_class$gene[match(mnar_idx, miss_class$feature_id)],
  feature_id = mnar_idx,
  pre_mean = rowMeans(mat[mnar_idx, , drop = FALSE], na.rm = TRUE),
  post_mean = rowMeans(mat_imp[mnar_idx, , drop = FALSE]),
  pre_sd = apply(mat[mnar_idx, , drop = FALSE], 1, sd, na.rm = TRUE),
  pct_miss = miss_class$pct_miss[match(mnar_idx, miss_class$feature_id)],
  imputation_reliable = miss_class$imputation_reliable[match(mnar_idx, miss_class$feature_id)]
) |>
  mutate(
    shift = post_mean - pre_mean,
    effect_d = shift / pre_sd
  )

imputed_df <- bind_cols(
  ann |>
    select(uniprot_id, protein, gene, description),
  as_tibble(mat_imp)
)
write_csv(imputed_df, file.path(cfg$data_dir, "01_imputed.csv"))

mask_df <- bind_cols(
  ann |>
    select(uniprot_id, gene, feature_id),
  as_tibble(was_na)
)
write_csv(mask_df, file.path(cfg$data_dir, "07_imputation_mask.csv"))
write_csv(miss_class, file.path(cfg$data_dir, "02_mar_mnar_classification.csv"))
write_csv(mnar_audit, file.path(cfg$data_dir, "08_mnar_imputation_audit.csv"))

dal <- dal_norm
mat_imp_uid <- mat_imp
rownames(mat_imp_uid) <- ann$uniprot_id
dal$data <- mat_imp_uid
dal$annotation <- dal$annotation |>
  left_join(
    miss_class |>
      select(uniprot_id, n_miss, pct_miss,
        miss_classification = classification,
        imputation_reliable),
    by = "uniprot_id"
  )
# Re-align $annotation rows to $data row order. dal_norm$annotation comes in
# normalization order; mat_imp_uid was reordered by gene_order for missForest
# determinism. Without this match() step the saved DAList has same set of
# proteins in $data and $annotation but different row positions.
dal$annotation <- dal$annotation[
  match(rownames(dal$data), dal$annotation$uniprot_id), , drop = FALSE]
rownames(dal$annotation) <- dal$annotation$uniprot_id
stopifnot(identical(rownames(dal$data), dal$annotation$uniprot_id))
saveRDS(dal, file.path(cfg$data_dir, "01_DAList_imputed.rds"))

summary_lines <- c(
  sprintf("n_proteins = %d", nrow(mat)),
  sprintf("n_samples = %d", ncol(mat)),
  sprintf("pct_missing = %.2f", pct_miss),
  sprintf("n_complete = %d", n_complete),
  sprintf("n_mar_proteins = %d", n_mar),
  sprintf("n_mnar_proteins = %d", n_mnar),
  sprintf("n_mar_values = %d", mar_vals),
  sprintf("n_mnar_values = %d", mnar_vals),
  sprintf("pct_mar_values = %.1f", pct_mar_vals),
  "classification_method = 3-method consensus (kmeans + logistic + left-tail)",
  sprintf("n_unreliable = %d", sum(!miss_class$imputation_reliable)),
  "best_method = missForest",
  sprintf("oob_error = %.4f", oob)
)
writeLines(summary_lines, file.path(cfg$data_dir, "09_imputation_summary.txt"))

saveRDS(
  list(
    mat = mat,
    mat_imp = mat_imp,
    was_na = was_na,
    ann = ann,
    meta = meta,
    miss_class = miss_class,
    miss_by_group = miss_by_group,
    prot_pct = prot_pct,
    pct_miss = pct_miss,
    mnar_genes = mnar_idx,
    mnar_audit = mnar_audit,
    best = "missForest",
    oob_error = oob,
    n_mar_prots = n_mar,
    n_mnar_prots = n_mnar,
    mar_miss_vals = mar_vals,
    mnar_miss_vals = mnar_vals,
    total_miss_vals = total_vals,
    classification_method = "3-method consensus (kmeans + logistic + left-tail)",
    PAL_GT = PAL_GT,
    PAL_MAR = PAL_MAR,
    PAL_CLASS = PAL_CLASS
  ),
  file.path(cfg$data_dir, "00_report_intermediates.rds")
)

writeLines(capture.output(sessionInfo()), file.path(cfg$data_dir, "sessionInfo.txt"))

cat(sprintf(
  "Done: missForest | %d proteins x %d samples | OOB = %.4f\n",
  nrow(mat_imp), ncol(mat_imp), oob
))
