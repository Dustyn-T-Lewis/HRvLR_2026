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
  corrected = "01_Trunk/03_Quantify/c_data/DAList_corrected.rds",
  limpa = "01_Trunk/03_Quantify/c_data/DAList_limpa.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
# Rank tests and correlations take no covariate, so they read limpa's complete blood-corrected
# matrix; observed_share marks how much of each protein was measured rather than estimated.
dal <- readRDS(paths[["corrected"]])
# Every protein carries one distinct gene symbol.
feature_label <- with(dal$annotation, set_names(gene, uniprot_id))
stopifnot(!anyNA(feature_label), !anyDuplicated(feature_label))
# Display order for the overviews: average-linkage clustering on 1 - Pearson r of limpa's matrix.
complete <- readRDS(paths[["limpa"]])$data

associate_step(
  stage = here("04_Continuous_no_group", "01_Associate", "01_Proteins"),
  inputs = inputs, paths = paths,
  values = dal$data, observed_share = rowMeans(dal$observed), meta = as_tibble(dal$metadata),
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
