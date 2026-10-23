library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "protein_grid.R"))

inputs <- c(contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx")
results <- read_sheet(inputs[["contrasts"]], "protein_by_contrast")
summary <- read_sheet(inputs[["contrasts"]], "contrast_summary")
n_samples <- nrow(read_sheet(inputs[["contrasts"]], "sample_weights"))
contrasts <- unique(results$contrast)

per_row <- 12
holds <- paste(
  "Panel b: per protein and contrast, the exploratory HR - LR candidates drawn (Π < 0.05 in an",
  sprintf("HR - LR or interaction contrast; top %d per block by p)", per_row),
  "with biopsies observed, in every contrast."
)
candidates <- hr_lr_candidates(results)
panel_data <- protein_panel_data(results, candidates, per_row)

# How many proteins each of these contrasts calls at p < 0.05, against what chance alone gives.
calls <- filter(summary, contrast_panel(contrast) %in% c("Between groups", "Interaction"))
chance <- calls$n_tested[1] * alpha_nominal
chance_range <- paste(range(calls$n_nominal), collapse = " to ")
fdr_text <- if (sum(calls$n_fdr_05) == 0) "none" else sum(calls$n_fdr_05)

panel <- panel_data |>
  protein_grid(
    contrasts,
    title = "Exploratory HR - LR candidates",
    subtitle = sprintf(
      paste0(
        "Π < %s in an HR - LR or interaction contrast; top %d per block by p\n",
        "p < %s: %s proteins per contrast, %s by chance; %s at FDR < %s\n",
        "(n/%d): biopsies where the protein was observed"
      ),
      alpha_nominal, per_row, alpha_nominal, chance_range, round(chance), fdr_text, alpha_fdr,
      n_samples
    ),
    n_samples = n_samples, limit = protein_colour_limit(results)
  )

save_panel(panel, "F2", "b_hr_lr_candidates", 95, 110)
