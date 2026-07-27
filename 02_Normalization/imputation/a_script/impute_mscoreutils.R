#!/usr/bin/env Rscript
# Imputation arm MsCoreUtils hybrid (exploratory; never the primary DEP input).
#
# Mechanism-aware hybrid, each tool used as designed: imputeLCMD::model.Selector()
# classifies each protein MAR (1) vs MNAR (0) (Lazar 2016), then
# MsCoreUtils::impute_matrix(method = "mixed") routes the MAR subset through kNN
# and the MNAR subset through QRILC (left-censored).

pacman::p_load(here, MsCoreUtils, imputeLCMD)
set.seed(42)

data_dir <- here("02_Normalization", "imputation", "c_data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cat(sprintf(
  "[mscoreutils] %d x %d | %.1f%% missing\n",
  nrow(mat), ncol(mat), mean(is.na(mat)) * 100
))

ms <- model.Selector(mat)
randna <- as.logical(ms[[1]]) # TRUE = MAR feature
cat(sprintf(
  "[mscoreutils] model.Selector split: %d MAR / %d MNAR features\n",
  sum(randna), sum(!randna)
))

# MARGIN is stated rather than left to impute_mixed's c(1, 1) default, because
# the second element is a real choice here. QRILC's own default is MARGIN = 2,
# one truncated normal per sample, which is the per-sample limit-of-detection
# model it was designed for. That model is unidentifiable here: model.Selector
# calls only a handful of proteins MNAR, so a per-sample fit draws on almost no
# observations and several samples have none at all. MARGIN = 2 fails outright
# (qnorm returns NaN, then lm.fit errors). MARGIN = 1 fits per protein across
# samples instead, which is not QRILC's intended model but is the only
# estimable one at this MNAR count. The printout above carries the live split.
imp <- impute_matrix(mat,
  method = "mixed", randna = randna,
  mar = "knn", mnar = "QRILC", MARGIN = c(1, 1)
)
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(
  method = "MsCoreUtils mixed (imputeLCMD model.Selector)",
  mar = "knn", mnar = "QRILC",
  n_mar = sum(randna), n_mnar = sum(!randna)
)
saveRDS(dal, file.path(data_dir, "DAList_imputed_mscoreutils.rds"))
cat("[mscoreutils] done: hybrid (knn/QRILC) -> DAList_imputed_mscoreutils.rds\n")
