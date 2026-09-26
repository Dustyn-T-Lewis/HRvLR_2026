# missForest copy of the normalised matrix, for the methods that need no missing values: fry,
# singscore and WGCNA. The limma fits read the unimputed matrix.

suppressPackageStartupMessages({
  library(here)
  library(missForest)
  library(tibble)
  library(writexl)
})

inputs <- c(normalized = "01_Preprocess/02_Normalization/c_data/DAList_normalized.rds")
paths <- vapply(inputs, here, character(1))
stopifnot(file.exists(paths))
dal <- readRDS(paths[["normalized"]])
out <- here("01_Preprocess", "03_Imputation", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
mat <- as.matrix(dal$data)

params <- list(method = "missForest", backend = "ranger", maxiter = 10, ntree = 100)
# Rows are sorted before the seeded fit so the result reproduces; ranger is named so a change in
# missForest's default cannot swap engines.
set.seed(42)
ord <- order(rownames(mat))
mf <- missForest(
  t(mat[ord, ]),
  maxiter = params$maxiter, ntree = params$ntree, backend = params$backend, verbose = FALSE
)
imputed <- t(mf$ximp)[rownames(mat), ]
dimnames(imputed) <- dimnames(mat)
stopifnot(!anyNA(imputed), identical(dim(imputed), dim(mat)))

dal$data <- imputed
dal$imputation <- c(params, oob_error = unname(mf$OOBerror[1]))
saveRDS(dal, file.path(out, "DAList_imputed.rds"), compress = "xz")
sheets <- list(
  imputation = as_tibble(c(
    params,
    pct_imputed = 100 * mean(is.na(mat)), oob_nrmse = unname(mf$OOBerror[1])
  )),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    dplyr::select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "missForest settings, share of values imputed, out-of-bag error.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "03_impute.xlsx"))
