# F5 a: every module's eigengene in the nine contrasts and its links to fibre area, muscle area and
# leg press, drawn by style/score_heatmap.R. Rows run largest module first, named by leading set.
library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "module_labels.R"))
source(here("04_Figures", "style", "score_heatmap.R"))

inputs <- c(
  module_tests = "02_HRvLR/03_Module_Contrasts/c_data/module_contrasts.rds",
  group_trait_slopes = "02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds",
  characterisation = characterisation_input
)
labels <- read_module_labels()
rows <- tibble(
  feature = labels$module, label = unname(module_identity(labels, width = 40)),
  n = labels$n_proteins
)
cells <- heatmap_cells(
  rename(readRDS(here(inputs[["module_tests"]]))$eigengene_tests, feature = module),
  readRDS(here(inputs[["group_trait_slopes"]])), "modules"
)
stopifnot(setequal(cells$feature, rows$feature))
panel_data <- cells |>
  inner_join(rename(rows, module = feature), by = c(feature = "module")) |>
  transmute(module = feature, label, n_proteins = n, section, column, key, scale, value, p, fdr)

panel <- score_heatmap(
  cells, rows, "proteins",
  title = "Modules: change, group difference and links to growth",
  subtitle = paste(
    "Green-purple: limma r-equivalent (contrasts; HR - LR slope); blue-red: Spearman rho within",
    "each group of 8. Outline: solid FDR < 0.05, dashed p < 0.05"
  )
)
save_panel(panel, "F5", "a_heatmap", 250, 95)
