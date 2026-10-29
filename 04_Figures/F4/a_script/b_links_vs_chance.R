# How many features track each trait: the share of proteins, set scores and module eigengenes at
# p < 0.05 in the change-score screen of 03_Continuous/01_Trait_Associations, per window and trait,
# against the 5% chance gives, with the count at FDR < 0.05.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "share_grid.R"))

inputs <- c(associations = "03_Continuous/01_Trait_Associations/c_data/trait_associations.rds")
associations <- readRDS(here(inputs[["associations"]]))

panel_data <- associations$counts |>
  filter(analysis == "change score") |>
  separate_wider_delim(screen, ", ", names = c("window", "trait")) |>
  mutate(n_expected = alpha_nominal * n_tested) |>
  select(feature_type, window, trait, n_tested, n_nominal, n_expected, n_fdr_05, share_nominal)

check_counts(panel_data, associations$change_score)

panel <- share_grid(
  panel_data,
  title = "How many features track each trait",
  subtitle = paste0(
    "Number: % of features at p < 0.05 (white: 5, chance); bold: count at FDR < 0.05\n",
    "Level windows pair a time point's values with trait levels; the rest with % change"
  )
)

save_panel(panel, "F4", "b_links_vs_chance", 165, 70)
