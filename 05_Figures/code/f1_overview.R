# Figure 1 for one grouping: who is compared (panel a), then the whole-proteome picture before any
# single result. It reads saved outputs only. The one computation is the phenotype PCA, a projection
# with no test; every p value drawn comes from steps 02 and 03 of 03_Group_comparisons.
#
# panel a     HRvLR: the BioRender design panel. Clusterings: each person in the space they were
#             clustered in, with the gap statistic's K, the K used and each group's bootstrap
#             Jaccard
# panel b     proteins per contrast and direction at Π < 0.05, p < 0.05 and BH < 0.05
# panel c     proteome PCA and phenotype PCA, with PERMANOVA
# panel d     selected gene sets across the nine contrasts (fgsea NES and FDR)
# panels e, f each group's training and acute sets, ranked in the other group (fry)
#
# build_group_overview() (00_Setup/Shared_functions/group_overview_figure.R) assembles panels b to
# f, checks every plotted number against its source table, and writes
# 05_Figures/<grouping>/F1_Overview/b_reports/f1_overview.pdf and c_data/f1_overview.xlsx.

library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_functions", "figure_style.R"))
source(here("00_Setup", "Shared_functions", "figure_panels.R"))
source(here("00_Setup", "Shared_functions", "group_overview_figure.R"))
source(here("00_Setup", "Shared_functions", "proteome_inputs.R"))

# What each clustering was built on, for panel a's subtitle.
clustered_on <- c(
  Muscle_size = "the % change in fibre area (mixed, type I, type II) and mCSA",
  Muscle_strength = "the % change in leg-press and leg-extension 1RM",
  Proteome_T1 = "the T1 module eigengenes", Proteome_T2 = "the T2 module eigengenes",
  Proteome_T3 = "the T3 module eigengenes",
  Proteome_training = "the training change in module eigengenes",
  Proteome_acute = "the acute change in module eigengenes"
)

f1_overview <- function(grouping) {
  comparisons <- file.path("03_Group_comparisons", grouping, "c_data")
  inputs <- c(
    contrasts = file.path(comparisons, "02_protein_contrasts.xlsx"),
    set_tests = file.path(comparisons, "03_set_tests.rds"),
    gene_sets = "01_Processing/04_Gene_Sets/c_data/gene_sets.rds",
    phenotype = "00_Setup/Input/phenotype.csv",
    if (grouping == "HRvLR") {
      c(design = "05_Figures/HRvLR/F1_Overview/assets/design.png")
    } else {
      c(
        clustering = file.path("02_Groups", grouping, "c_data/clustering.xlsx"),
        design = file.path(comparisons, "01_design.rds")
      )
    }
  )
  paths <- map_chr(inputs, here)
  stopifnot(file.exists(paths))
  out_dir <- here("05_Figures", grouping, "F1_Overview")
  if (grouping == "HRvLR") {
    # Panel a is the BioRender design panel, pinned under its title so it sits level with b and c.
    panel_design <- ggplot() +
      annotation_custom(grid::rasterGrob(
        png::readPNG(paths[["design"]]),
        y = grid::unit(1, "npc"), just = "top", interpolate = TRUE
      )) +
      labs(
        title = "Design and contrasts",
        subtitle = "16 adults, 45 biopsies; contrast colours carry through b to f"
      ) +
      theme_void(base_size = 7) +
      theme(
        plot.title = figure_theme$plot.title, plot.subtitle = figure_theme$plot.subtitle,
        plot.title.position = "plot"
      )
    return(build_group_overview(
      panel_design, inputs, paths,
      groups = c("HR", "LR"), out_dir = out_dir, stem = "f1_overview",
      note = "Panel a is the BioRender design panel (assets/design.png); it reads no data."
    ))
  }
  read_clustering <- \(sheet) readxl::read_excel(paths[["clustering"]], sheet)
  is_proteome <- startsWith(grouping, "Proteome_")
  groups <- if (is_proteome) c("C1", "C2") else c("Hi", "Lo")
  space <- read_clustering("space")
  if (is_proteome) space <- mutate(space, group = paste0("C", cluster))
  clustering <- panel_clustering(
    space, read_clustering("space_variance"), read_clustering("choice"),
    read_clustering("stability"), groups,
    clustered_on = clustered_on[[grouping]]
  )
  panel_a <- clustering$plot
  validity <- NULL
  if (is_proteome) {
    # Each contrast's label carries its validity against the clustering: circ. (circular, built
    # from the same time points) or caut. (caution, shares a time point).
    input <- sub("Proteome_", "", grouping)
    panel_a <- panel_a + labs(title = sprintf("Who is in each group, clustered at %s", input))
    contrasts <- as.character(overview_contrasts(groups)$contrast)
    validity <- set_names(validity(input, contrasts), contrasts)
  }
  build_group_overview(
    panel_a, inputs, paths,
    groups = groups, out_dir = out_dir, stem = "f1_overview",
    a_sheets = list(clustering = clustering$data),
    a_holds = if (is_proteome) {
      "Panel a: each person's cluster, HR/LR and PC1 and PC2 of the eigengenes."
    } else {
      "Panel a: each person's cluster, HR/LR and PC1 and PC2 of the scaled traits."
    },
    validity = validity,
    group_of = readRDS(paths[["design"]])$group_of
  )
}
