library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "protein_grid.R"))

inputs <- c(contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx")
results <- read_sheet(inputs[["contrasts"]], "protein_by_contrast")
n_samples <- nrow(read_sheet(inputs[["contrasts"]], "sample_weights"))
contrasts <- unique(results$contrast)

# A responder has a Π hit in a group's training or acute contrast. Each group's evidence is the
# smallest p among its hits; shared proteins rank on the weaker group's.
per_row <- 10
holds <- paste(
  "Panel a: per protein and contrast, the responders drawn (Π < 0.05 in a group's training or",
  sprintf("acute contrast; top %d per block by p)", per_row),
  "with biopsies observed, logFC, 95% CI, p, FDR and Π in every contrast."
)
within <- c("Training_HR", "Acute_HR", "Training_LR", "Acute_LR")
block_levels <- c("Shared", "HR only", "LR only")
chosen <- results |>
  filter(contrast %in% within, pi_score < alpha_nominal) |>
  mutate(group = sub(".*_", "", contrast)) |>
  summarise(best = min(p), .by = c(uniprot_id, group)) |>
  pivot_wider(names_from = group, values_from = best) |>
  mutate(
    block = factor(case_when(
      !is.na(HR) & !is.na(LR) ~ block_levels[1],
      !is.na(HR) ~ block_levels[2],
      .default = block_levels[3]
    ), block_levels),
    selection_p = pmax(HR, LR, na.rm = TRUE)
  )
panel_data <- protein_panel_data(results, chosen, per_row)

panel <- panel_data |>
  protein_grid(
    contrasts,
    title = "Proteins that respond to training or the bout",
    subtitle = sprintf(
      paste0(
        "Π < %s in a group's training or acute contrast; top %d per block by p\n",
        "(shared: the weaker group's p); (n/%d): biopsies where the protein was observed"
      ),
      alpha_nominal, per_row, n_samples
    ),
    n_samples = n_samples, limit = protein_colour_limit(results)
  )

save_panel(panel, "F2", "a_responders", 95, 130)
