#!/usr/bin/env Rscript
# Imputation arm C: missForest (exploratory; never the primary DEP input).
#
# Random-forest imputation (Stekhoven & Buhlmann 2012), mechanism-agnostic: one
# model for all missing values. Rows are sorted before the stochastic fit so the
# result is reproducible. This is the arm downstream QC/figures read by default.

pacman::p_load(proteoDA, here, missForest)
set.seed(42)

data_dir <- here("02_Normalization", "imputation", "c_data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cat(sprintf(
  "[missforest] %d x %d | %.1f%% missing\n",
  nrow(mat), ncol(mat), mean(is.na(mat)) * 100
))

ord <- order(rownames(mat)) # deterministic order before the stochastic fit
mf <- missForest(t(mat[ord, ]), maxiter = 10, ntree = 100, verbose = TRUE)
imp <- t(mf$ximp)[rownames(mat), ]
dimnames(imp) <- dimnames(mat)
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))
cat(sprintf("[missforest] OOB error: %.4f\n", mf$OOBerror[1]))

dal$data <- imp
dal$imputation <- list(
  method = "missForest", maxiter = 10, ntree = 100,
  oob_error = as.numeric(mf$OOBerror[1])
)
saveRDS(dal, file.path(data_dir, "DAList_imputed_missforest.rds"))
cat("[missforest] done -> DAList_imputed_missforest.rds\n")
