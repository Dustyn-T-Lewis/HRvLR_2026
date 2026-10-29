# S3 e: the qualifying Hallmark sets as singscore set scores, in the F5 a layout
# (style/score_heatmap.R). Hallmark is chosen in advance, so the rows are not picked by the tests
# drawn; they are ordered by clustering their 36 values.
library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "score_heatmap.R"))

inputs <- c(
  set_tests = "02_HRvLR/02_Set_Contrasts/c_data/set_contrasts.rds",
  group_trait_slopes = "02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds",
  gene_sets = "01_Processing/04_Gene_Sets/c_data/gene_sets.rds"
)
hallmark <- readRDS(here(inputs[["gene_sets"]]))$set_catalog |>
  filter(collection == "Hallmark", qualifies)
cells <- heatmap_cells(
  rename(readRDS(here(inputs[["set_tests"]]))$score_tests, feature = set_id),
  readRDS(here(inputs[["group_trait_slopes"]])), "sets"
) |>
  filter(feature %in% hallmark$set_id)
profile <- cells |>
  select(feature, key, value) |>
  tidyr::pivot_wider(names_from = key, values_from = value) |>
  tibble::column_to_rownames("feature") |>
  as.matrix()
rows <- tibble(feature = rownames(profile)[hclust(dist(profile))$order]) |>
  left_join(select(hallmark, feature = set_id, pathway, n = measured_size), by = "feature") |>
  mutate(label = pathway_label(pathway, 40))
stopifnot(nrow(rows) == nrow(hallmark), setequal(cells$feature, rows$feature))
panel_data <- cells |>
  inner_join(select(rows, feature, label, n), by = "feature") |>
  transmute(set_id = feature, label, measured_size = n, section, column, key, scale, value, p, fdr)

panel <- score_heatmap(
  cells, rows, "measured genes",
  title = "Hallmark set scores: change, group difference and links to growth",
  subtitle = sprintf(
    paste(
      "The %d Hallmark sets with 15 to 500 measured genes, as singscore scores; columns and",
      "colours as Figure 5a; rows clustered by their values"
    ),
    nrow(rows)
  )
)
panel_height <- 4.2 * nrow(rows) + 45
save_panel(panel, "S3", "e_hallmark_heatmap", 250, panel_height)
