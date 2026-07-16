#!/usr/bin/env Rscript
# Imputation arm MsCoreUtils hybrid (exploratory; never the primary DEP input).
#
# Mechanism-aware hybrid, each tool used as designed: imputeLCMD::model.Selector()
# classifies each protein MAR (1) vs MNAR (0) (Lazar 2016), then
# MsCoreUtils::impute_matrix(method = "mixed") routes the MAR subset through kNN
# and the MNAR subset through QRILC (left-censored).

pacman::p_load(proteoDA, here, MsCoreUtils, imputeLCMD)
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

imp <- impute_matrix(mat,
  method = "mixed", randna = randna,
  mar = "knn", mnar = "QRILC"
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
