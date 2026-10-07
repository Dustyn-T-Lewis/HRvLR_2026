# Run once per proteome clustering (T1, T2, T3, training change, acute change) on that
# clustering's groups C1 and C2: a copy of 02_HRvLR_groups/05_Classify/01_Proteins. Each input
# writes to its own folder; stack_inputs() then joins them into c_data/all_inputs.xlsx with a
# validity tag. How well each protein separates the groups, on the five between-group tasks of
# classify_step().

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
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

run_input <- function(input) {
  inputs <- c(
    corrected = "03_Cluster_groups/Proteome/01_Cluster/c_data/inputs/{input}/DAList_corrected.rds",
    limpa = "03_Cluster_groups/Proteome/01_Cluster/c_data/inputs/{input}/DAList_limpa.rds"
  )
  inputs <- sub("{input}", input, inputs, fixed = TRUE)
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

  classify_step(
    stage = here("03_Cluster_groups", "Proteome", "06_Classify", "01_Proteins", input),
    inputs = inputs, paths = paths,
    values = dal$data, observed_share = rowMeans(dal$observed), meta = dal$metadata,
    feature_label = feature_label, feature_noun = "proteins", feature_unit = "protein",
    value_label = "log2 abundance",
    row_groups = list(
      proteins = rownames(complete)[hclust(as.dist(1 - cor(t(complete))), "average")$order]
    ),
    row_order = paste(
      "Rows are ordered by average-linkage clustering (1 - Pearson r) of limpa's",
      "matrix, for display only."
    ),
    groups = c("C1", "C2")
  )
}
purrr::walk(inputs_used, run_input)
stack_inputs(here("03_Cluster_groups", "Proteome", "06_Classify", "01_Proteins"))
