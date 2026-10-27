# How many features link to a trait differently in HR and LR: the share of proteins, set scores and
# module eigengenes whose HR slope minus LR slope (02_HRvLR/07_Group_Trait_Slopes) sits at
# p < 0.05, per window and trait, against the 5% chance gives, with the count at FDR < 0.05.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "share_grid.R"))

inputs <- c(group_trait_slopes = "02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds")
slopes <- readRDS(here(inputs[["group_trait_slopes"]]))
holds <- paste(
  "Panel a: features whose HR slope minus LR slope sits at p < 0.05 and FDR < 0.05, per feature",
  "type, window and trait, against the 5% chance gives."
)

panel_data <- slopes$counts |>
  filter(test == "group_difference") |>
  mutate(n_expected = alpha_nominal * n_tested) |>
  select(feature_type, window, trait, n_tested, n_nominal, n_expected, n_fdr_05, share_nominal)

check_counts(panel_data, slopes$group_difference)
stopifnot(!"comp_hypertrophy" %in% panel_data$trait)

panel <- share_grid(
  panel_data,
  title = "Does a feature follow a trait differently in HR and LR?",
  subtitle = paste0(
    "HR slope minus LR slope (limma, feature ~ group x trait). Number: % of features at ",
    "p < 0.05 (white: 5, chance);\nbold below: count at FDR < 0.05. The composite score defined ",
    "the groups and is left out"
  )
)

save_panel(panel, "S4", "a_slope_counts", 250, 65)
