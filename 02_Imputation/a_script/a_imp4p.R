#!/usr/bin/env Rscript
# Imputation arm A: imp4p (exploratory; never the primary DEP input).
#
# imp4p (Giai Gianetto et al.) is built for label-free proteomics with a mixture
# of MCAR/MNAR missingness and estimates the mechanism itself, so we pass only
# the matrix and the Group_Time factor. Output feeds QC and sensitivity work; the
# primary DEP stays on the non-imputed normalized matrix.

library(proteoDA)
library(imp4p)

set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

norm_rds <- "01_normalization/c_data/03_DAList_normalized.rds"
data_dir <- "02_Imputation/c_data"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(norm_rds)
mat <- as.matrix(dal$data)
cond <- factor(dal$metadata$Group_Time[match(colnames(mat), dal$metadata$Col_ID)])
cat(sprintf(
  "[imp4p] %d x %d | %.1f%% missing | conditions: %s\n",
  nrow(mat), ncol(mat), mean(is.na(mat)) * 100,
  paste(levels(cond), collapse = "/")
))

# mle for the MCAR part, igcda for the MNAR part; imp4p decides the per-protein mix
imp <- as.matrix(impute.mi(
  tab = mat, conditions = cond,
  methodMCAR = "mle", methodMNAR = "igcda",
  progress.bar = FALSE
))
dimnames(imp) <- dimnames(mat)
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(
  method = "imp4p::impute.mi", methodMCAR = "mle",
  methodMNAR = "igcda"
)
saveRDS(dal, file.path(data_dir, "DAList_imputed_imp4p.rds"))
cat(sprintf(
  "[imp4p] done: imputed %d cells -> DAList_imputed_imp4p.rds\n",
  sum(is.na(mat))
))
