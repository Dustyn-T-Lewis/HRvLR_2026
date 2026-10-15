library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "protein_grid.R"))

inputs <- c(contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx")
results <- readxl::read_excel(here(inputs[["contrasts"]]), "protein_by_contrast")
n_samples <- nrow(readxl::read_excel(here(inputs[["contrasts"]]), "sample_weights"))
contrasts <- unique(results$contrast)
stopifnot(length(contrasts) == 9)

# A responder is a protein at Π < 0.05 in a group's training or acute contrast. Each group's
# evidence is its smallest Π over those two; shared proteins rank on the weaker group's.
per_row <- 10
holds <- paste(
  "Panel a: per protein and contrast, the responders drawn (Π < 0.05 in a group's training or",
  sprintf("acute contrast; top %d per block by Π)", per_row),
  "with biopsies observed, logFC, 95% CI, p, FDR and Π in all nine contrasts."
)
within <- c("Training_HR", "Acute_HR", "Training_LR", "Acute_LR")
block_levels <- c("Shared", "HR only", "LR only")
chosen <- results |>
  filter(contrast %in% within) |>
  mutate(group = sub(".*_", "", contrast)) |>
  summarise(best = min(pi_score), .by = c(uniprot_id, group)) |>
  pivot_wider(names_from = group, values_from = best) |>
  mutate(
    block = case_when(
      HR < alpha_nominal & LR < alpha_nominal ~ block_levels[1],
      HR < alpha_nominal ~ block_levels[2],
      LR < alpha_nominal ~ block_levels[3]
    ),
    selection_pi = case_when(
      block == block_levels[1] ~ pmax(HR, LR), block == block_levels[2] ~ HR,
      .default = LR
    )
  ) |>
  filter(!is.na(block)) |>
  slice_min(selection_pi, n = per_row, by = block, with_ties = FALSE) |>
  mutate(block = factor(block, block_levels), row_order = rank(selection_pi), .by = block)

panel_data <- results |>
  inner_join(select(chosen, uniprot_id, block, row_order, selection_pi), by = "uniprot_id") |>
  select(
    block, row_order, uniprot_id, gene, n_observed, contrast, logFC, logFC_lower, logFC_upper,
    p, fdr, pi_score, selection_pi
  ) |>
  arrange(block, row_order, match(contrast, contrasts))

check <- inner_join(panel_data, results, by = c("uniprot_id", "contrast"))
stopifnot(
  nrow(panel_data) == 9 * nrow(chosen), nrow(check) == nrow(panel_data),
  all(check$logFC.x == check$logFC.y), all(check$p.x == check$p.y),
  all(check$fdr.x == check$fdr.y), all(check$pi_score.x == check$pi_score.y),
  all(chosen$selection_pi < alpha_nominal), max(count(chosen, block)$n) <= per_row
)

panel <- panel_data |>
  rename(order = row_order) |>
  protein_grid(
    contrasts,
    title = "Proteins that respond to training or the bout",
    subtitle = sprintf(
      paste0(
        "Π < %s in a group's training or acute contrast\n",
        "top %d per block by Π (shared: the weaker group's)\n",
        "(n/%d): biopsies where the protein was observed"
      ),
      alpha_nominal, per_row, n_samples
    ),
    n_samples = n_samples
  )

ggsave(
  here("04_Figures", "F2", "b_reports", "panels", "a_responders.pdf"), panel,
  width = 95, height = 130, units = "mm", device = grDevices::cairo_pdf
)
