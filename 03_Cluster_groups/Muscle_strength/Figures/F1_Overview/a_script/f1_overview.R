# Figure 1 for the strength grouping: who the clusters are, then the whole-proteome picture
# with the clusters in place of HR and LR. Reads saved outputs only; the one computation is the
# phenotype PCA, a projection with no test.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source(here("R", "helpers.R"))
source(here("R", "figure_theme.R"))
source(here("R", "figure_panels.R"))
source(here("R", "figure_overview.R"))

folder <- "03_Cluster_groups/Muscle_strength"
inputs <- c(
  contrasts = file.path(folder, "03_Protein_Contrasts/c_data/protein_contrasts.xlsx"),
  set_tests = file.path(folder, "04_Set_Contrasts/c_data/set_tests.rds"),
  gene_sets = "01_Processing/04_Gene_Sets/c_data/gene_sets.rds",
  phenotype = "00_Input/phenotype.csv",
  clustering = file.path(folder, "01_Clustering/c_data/clustering.xlsx"),
  design = file.path(folder, "02_Design_Contrasts/c_data/design.rds")
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
read_clustering <- \(sheet) readxl::read_excel(paths[["clustering"]], sheet)
groups <- c("Hi", "Lo")

clustering <- panel_clustering(
  read_clustering("space"), read_clustering("space_variance"), read_clustering("choice"),
  read_clustering("stability"), groups,
  clustered_on = "the % change in leg-press and leg-extension 1RM"
)
figure_overview(
  clustering$plot, inputs, paths,
  groups = groups, out_dir = here(folder, "Figures", "F1_Overview"), stem = "f1_overview",
  a_sheets = list(clustering = clustering$data),
  a_holds = "Panel a: each person's cluster, HR/LR and PC1 and PC2 of the scaled traits.",
  group_of = readRDS(paths[["design"]])$group_of
)
message("wrote f1_overview.pdf and f1_overview.xlsx")
