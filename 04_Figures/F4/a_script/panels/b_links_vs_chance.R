# How many features track each trait: the share of proteins, set scores and module eigengenes at
# p < 0.05 in the change-score screen of 03_Continuous/01_Trait_Associations, per window and trait,
# against the 5% chance gives, with the count at FDR < 0.05.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "share_grid.R"))

inputs <- c(associations = "03_Continuous/01_Trait_Associations/c_data/associations.rds")
associations <- readRDS(here(inputs[["associations"]]))
holds <- paste(
  "Panel b: features at p < 0.05 and FDR < 0.05 per feature type, window and trait (change score),",
  "against the 5% chance gives."
)

panel_data <- associations$counts |>
  filter(analysis == "change score") |>
  separate_wider_delim(screen, ", ", names = c("window", "trait")) |>
  mutate(n_expected = alpha_nominal * n_tested) |>
  select(feature_type, window, trait, n_tested, n_nominal, n_expected, n_fdr_05, share_nominal)

recount <- associations$change_score |>
  summarise(
    n_tested = sum(!is.na(p)), n_nominal = sum(p < alpha_nominal, na.rm = TRUE),
    n_fdr_05 = sum(fdr < alpha_fdr, na.rm = TRUE), .by = c(feature_type, window, trait)
  )
stopifnot(
  nrow(panel_data) == nrow(recount),
  identical(
    arrange(
      select(panel_data, feature_type, window, trait, n_tested, n_nominal, n_fdr_05),
      feature_type, window, trait
    ),
    arrange(recount, feature_type, window, trait)
  )
)

panel <- share_grid(
  panel_data,
  title = "How many features track each trait",
  subtitle = paste0(
    "Number: % of features at p < 0.05 (white: 5, chance); bold: count at FDR < 0.05\n",
    "Level windows pair a time point's values with trait levels; the rest with % change"
  )
)

ggsave(
  here("04_Figures", "F4", "b_reports", "panels", "b_links_vs_chance.pdf"), panel,
  width = 105, height = 70, units = "mm", device = cairo_pdf
)
