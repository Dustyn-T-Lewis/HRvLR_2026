library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "protein_grid.R"))

inputs <- c(contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx")
read_sheet <- \(sheet) readxl::read_excel(here(inputs[["contrasts"]]), sheet)
results <- read_sheet("protein_by_contrast")
summary <- read_sheet("contrast_summary")
n_samples <- nrow(read_sheet("sample_weights"))
contrasts <- unique(results$contrast)
stopifnot(length(contrasts) == 9)

# A candidate is a protein at Π < 0.05 in an HR - LR contrast (T1, T2, T3) or an interaction. It
# sits in the block of the contrast where its Π is smallest and ranks on that Π.
per_row <- 12
holds <- paste(
  "Panel b: per protein and contrast, the exploratory HR - LR candidates drawn (Π < 0.05 in an",
  sprintf("HR - LR or interaction contrast; top %d per block by Π)", per_row),
  "with biopsies observed, in all nine contrasts."
)
hr_lr_contrasts <- contrasts[contrast_panel(contrasts) %in% c("Between groups", "Interaction")]
chosen <- results |>
  filter(contrast %in% hr_lr_contrasts) |>
  slice_min(pi_score, n = 1, by = uniprot_id, with_ties = FALSE) |>
  filter(pi_score < alpha_nominal) |>
  transmute(uniprot_id,
    block = contrast_panel(contrast), selection_pi = pi_score,
    selected_in = contrast
  ) |>
  slice_min(selection_pi, n = per_row, by = block, with_ties = FALSE) |>
  mutate(block = droplevels(block), row_order = rank(selection_pi), .by = block)

panel_data <- results |>
  inner_join(
    select(chosen, uniprot_id, block, row_order, selected_in, selection_pi),
    by = "uniprot_id"
  ) |>
  select(
    block, row_order, uniprot_id, gene, n_observed, contrast, logFC, logFC_lower, logFC_upper,
    p, fdr, pi_score, selected_in, selection_pi
  ) |>
  arrange(block, row_order, match(contrast, contrasts))

check <- inner_join(panel_data, results, by = c("uniprot_id", "contrast"))
stopifnot(
  nrow(panel_data) == 9 * nrow(chosen), nrow(check) == nrow(panel_data),
  all(check$logFC.x == check$logFC.y), all(check$p.x == check$p.y),
  all(check$fdr.x == check$fdr.y), all(check$pi_score.x == check$pi_score.y),
  all(chosen$selection_pi < alpha_nominal), max(count(chosen, block)$n) <= per_row
)

# How many proteins each of these contrasts calls at p < 0.05, against what chance alone gives.
calls <- filter(summary, contrast %in% hr_lr_contrasts)
chance <- calls$n_tested[1] * alpha_nominal
chance_range <- paste(range(calls$n_nominal), collapse = " to ")
fdr_text <- if (sum(calls$n_fdr_05) == 0) "none" else sum(calls$n_fdr_05)

panel <- panel_data |>
  rename(order = row_order) |>
  protein_grid(
    contrasts,
    title = "Exploratory HR - LR candidates",
    subtitle = sprintf(
      paste0(
        "Π < %s in an HR - LR or interaction contrast; top %d per block by Π\n",
        "p < %s: %s proteins per contrast, %s by chance; %s at FDR < %s\n",
        "(n/%d): biopsies where the protein was observed. Seen in few,\n",
        "a protein's values are mostly modelled and its logFC can be large"
      ),
      alpha_nominal, per_row, alpha_nominal, chance_range, round(chance), fdr_text, alpha_fdr,
      n_samples
    ),
    n_samples = n_samples
  )

ggsave(
  here("04_Figures", "F2", "b_reports", "panels", "b_hr_lr_candidates.pdf"), panel,
  width = 95, height = 110, units = "mm", device = grDevices::cairo_pdf
)
