# Ties each gene set's singscore to phenotype by the sample model and the change scores in
# R/associate.R. Each collection is its own screen: BH runs within collection.

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
  singscore = "01_Trunk/05_Scores/c_data/singscore.rds",
  gene_sets = "01_Trunk/04_Gene_Sets/c_data/gene_sets.rds",
  proteins = "01_Trunk/03_Quantify/c_data/DAList_corrected.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
catalog <- readRDS(paths[["gene_sets"]])$set_catalog |> filter(qualifies)
values <- readRDS(paths[["singscore"]])$scores[catalog$set_id, ]
meta <- as_tibble(readRDS(paths[["proteins"]])$metadata)
stopifnot(identical(colnames(values), meta$sample_id))

associate_step(
  stage = here("04_Continuous_no_group", "01_Associate", "02_Sets"),
  inputs = inputs, paths = paths,
  values = values, meta = meta,
  feature_label = with(catalog, set_names(enrichVolcano::ev_clean_label(pathway), set_id)),
  feature_noun = "sets", feature_unit = "set", value_label = "singscore",
  # Display order within each collection: average-linkage clustering on 1 - Pearson r of
  # singscore.
  row_groups = split(catalog$set_id, factor(catalog$collection, unique(catalog$collection))) |>
    map(\(ids) ids[hclust(as.dist(1 - cor(t(values[ids, ]))), "average")$order]),
  row_order = paste(
    "Rows are ordered within collection by average-linkage clustering (1 - Pearson r)",
    "of singscore, for display only."
  ),
  collection_of = with(catalog, set_names(collection, set_id)),
  # ev_clean_label() breaks names every few words; panel headers re-wrap them at 50 characters.
  panel_name = \(label, collection) {
    paste0(collection, ": ", stringr::str_wrap(gsub("\n", " ", label), 50))
  }
)
