#!/usr/bin/env Rscript
# Imputation arm Perseus MNAR left-censoring (proteoDA::perseus_impute).
#
# Per-sample downshifted-Gaussian draw (Tyanova 2016), the MNAR-aware contrast to
# missForest's MAR. Exploratory only; never the reported DEP input. Comparison arm.

pacman::p_load(proteoDA, here)

data_dir <- here("02_Normalization", "imputation", "c_data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cat(sprintf(
  "[perseus] %d x %d | %.1f%% missing\n",
  nrow(mat), ncol(mat), mean(is.na(mat)) * 100
))

dal <- perseus_impute(dal, shift = 1.8, width = 0.3, robust = TRUE, seed = 42)
imp <- as.matrix(dal$data)
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$imputation <- list(method = "perseus", shift = 1.8, width = 0.3, robust = TRUE)
saveRDS(dal, file.path(data_dir, "DAList_imputed_perseus.rds"))
cat("[perseus] done -> DAList_imputed_perseus.rds\n")
