# Figure 1 for the proteome grouping, once per clustering input (T1, T2, T3, training change, acute
# change): who the clusters are, then the whole-proteome picture with C1 and C2 in place of HR and
# LR. Each contrast's label carries its validity against that clustering: circ. (circular, built
# from the same time points) or caut. (caution, shares a time point). Reads saved outputs only.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_functions", "figure_style.R"))
source(here("00_Setup", "Shared_functions", "figure_panels.R"))
source(here("00_Setup", "Shared_functions", "group_overview_figure.R"))
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

folder <- "03_Cluster_groups/Proteome"
groups <- c("C1", "C2")
clustered_on <- c(
  T1 = "the T1 module eigengenes", T2 = "the T2 module eigengenes",
  T3 = "the T3 module eigengenes", training = "the training change in module eigengenes",
  acute = "the acute change in module eigengenes"
)

figure_input <- function(input) {
  inputs <- c(
    contrasts = file.path(folder, "03_Protein_Contrasts", input, "c_data/protein_contrasts.xlsx"),
    set_tests = file.path(folder, "04_Set_Contrasts", input, "c_data/set_tests.rds"),
    gene_sets = "01_Processing/04_Gene_Sets/c_data/gene_sets.rds",
    phenotype = "00_Setup/Input/phenotype.csv",
    clustering = "02_Groups/Proteome/c_data/clustering.xlsx",
    design = file.path(folder, "02_Design_Contrasts", input, "c_data/design.rds")
  )
  paths <- map_chr(inputs, here)
  stopifnot(file.exists(paths))
  read_clustering <- \(sheet) {
    readxl::read_excel(paths[["clustering"]], sheet) |> filter(.data$input == .env$input)
  }
  space <- read_clustering("space") |> mutate(group = paste0("C", cluster))
  clustering <- panel_clustering(
    space, read_clustering("space_variance"), read_clustering("choice"),
    read_clustering("stability"), groups, clustered_on[[input]]
  )
  contrasts <- as.character(overview_contrasts(groups)$contrast)
  build_group_overview(
    clustering$plot + labs(title = sprintf("Who is in each group, clustered at %s", input)),
    inputs, paths,
    groups = groups, out_dir = here(folder, "Figures", "F1_Overview"),
    stem = paste0("f1_overview_", input),
    a_sheets = list(clustering = clustering$data),
    a_holds = "Panel a: each person's cluster, HR/LR and PC1 and PC2 of the eigengenes.",
    validity = set_names(validity(input, contrasts), contrasts),
    group_of = readRDS(paths[["design"]])$group_of
  )
}
walk(inputs_used, figure_input)
