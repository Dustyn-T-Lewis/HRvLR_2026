# Ties each protein to phenotype by the sample model and the change scores in R/associate.R.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(writexl)
  library(grid)
  library(ComplexHeatmap)
})
source(here("R", "helpers.R"))
source(here("R", "associate.R"))

inputs <- c(
  corrected = "01_Preprocess/03_Quantification/c_data/DAList_corrected.rds",
  limpa = "01_Preprocess/03_Quantification/c_data/DAList_limpa.rds",
  phenotype = "00_Input/phenotype.csv"
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

associate_step(
  stage = here("02_Differential_Expression", "04_Associate"), inputs = inputs, paths = paths,
  values = replace(dal$data, !dal$observed, NA), meta = as_tibble(dal$metadata),
  feature_label = feature_label, feature_noun = "proteins", feature_unit = "protein",
  value_label = "log2 abundance",
  row_groups = list(
    proteins = rownames(complete)[hclust(as.dist(1 - cor(t(complete))), "average")$order]
  ),
  row_order = paste(
    "Rows are ordered by average-linkage clustering (1 - Pearson r) of limpa's",
    "matrix, for display only."
  ),
  save_change_score = TRUE
)
