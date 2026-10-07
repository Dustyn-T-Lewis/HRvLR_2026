# Run once per proteome clustering (T1, T2, T3, training change, acute change) on that
# clustering's groups C1 and C2: a copy of 02_HRvLR_groups/06_Classification/02_Sets. Each input
# writes to its own folder; stack_inputs() then joins them into c_data/<content>_all_inputs.xlsx
# with a validity tag. How well each gene set's singscore separates the groups, on the five between-
# group tasks of classify_step(). Each collection is its own screen, so BH runs within task and
# collection.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ComplexHeatmap)
})
source(here("R", "helpers.R"))
source(here("R", "classify.R"))
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

run_input <- function(input) {
  inputs <- c(
    singscore = "03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/{input}/singscore.rds",
    gene_sets = "01_Processing/04_Gene_Sets/c_data/gene_sets.rds",
    proteins = "03_Cluster_groups/Proteome/01_Clustering/c_data/inputs/{input}/DAList_corrected.rds"
  )
  inputs <- sub("{input}", input, inputs, fixed = TRUE)
  paths <- map_chr(inputs, here)
  stopifnot(file.exists(paths))
  catalog <- readRDS(paths[["gene_sets"]])$set_catalog |> filter(qualifies)
  values <- readRDS(paths[["singscore"]])$scores[catalog$set_id, ]
  meta <- readRDS(paths[["proteins"]])$metadata
  stopifnot(identical(colnames(values), meta$sample_id))

  classify_step(
    stage = here("03_Cluster_groups", "Proteome", "06_Classification", "02_Sets", input),
    name = "set_classification", inputs = inputs, paths = paths,
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
    },
    groups = c("C1", "C2")
  )
}
purrr::walk(inputs_used, run_input)
stack_inputs(here("03_Cluster_groups", "Proteome", "06_Classification", "02_Sets"))
