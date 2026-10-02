# How well each protein separates the study's groups, on the eight tasks in R/classify.R.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ComplexHeatmap)
  library(writexl)
})
source(here("R", "helpers.R"))
source(here("R", "classify.R"))

inputs <- c(
  corrected = "01_Preprocess/03_Quantification/c_data/DAList_corrected.rds",
  limpa = "01_Preprocess/03_Quantification/c_data/DAList_limpa.rds"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
# Rank tests and correlations take no covariate, so they read the blood-corrected matrix, and
# only where a value was observed: limpa's estimates for missing values carry wide errors that
# these tests cannot weigh.
dal <- readRDS(paths[["corrected"]])
# Every protein carries one distinct gene symbol.
feature_label <- with(dal$annotation, set_names(gene, uniprot_id))
stopifnot(!anyNA(feature_label), !anyDuplicated(feature_label))
# Display order for the overviews: average-linkage clustering on 1 - Pearson r of limpa's matrix.
complete <- readRDS(paths[["limpa"]])$data

classify_step(
  stage = here("02_Differential_Expression", "03_Classify"), inputs = inputs, paths = paths,
  values = replace(dal$data, !dal$observed, NA), meta = dal$metadata,
  feature_label = feature_label, feature_noun = "proteins", feature_unit = "protein",
  value_label = "log2 abundance",
  row_groups = list(
    proteins = rownames(complete)[hclust(as.dist(1 - cor(t(complete))), "average")$order]
  ),
  row_order = paste(
    "Rows are ordered by average-linkage clustering (1 - Pearson r) of limpa's",
    "matrix, for display only."
  )
)
