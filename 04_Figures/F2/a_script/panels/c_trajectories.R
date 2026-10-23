library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "protein_grid.R"))

inputs <- c(
  contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx",
  corrected = "01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds"
)
results <- read_sheet(inputs[["contrasts"]], "protein_by_contrast")
corrected <- readRDS(here(inputs[["corrected"]]))
meta <- as_tibble(corrected$targets)

# The HR - LR candidates with the smallest p (the panel b rule without the per-block cap).
n_shown <- 4
holds <- paste(
  sprintf("Panel c: per biopsy, the blood-corrected value of the %d HR - LR candidates", n_shown),
  "with the smallest p, and whether a precursor was observed."
)
chosen <- hr_lr_candidates(results) |>
  slice_min(selection_p, n = n_shown, with_ties = FALSE) |>
  mutate(facet = sprintf(
    "%s (%d/%d)\n%s\nlogFC %.2f, Π = %s\n%s, %s",
    gene, as.integer(n_observed), ncol(corrected$E), contrast_label(selected_in), logFC,
    formatC(pi_score, format = "e", digits = 1), label_p(p), label_fdr(fdr)
  ))

panel_data <- corrected$E[chosen$uniprot_id, , drop = FALSE] |>
  as_tibble(rownames = "uniprot_id") |>
  tidyr::pivot_longer(-uniprot_id, names_to = "sample_id", values_to = "log2_abundance") |>
  inner_join(select(meta, sample_id, subject, group, timepoint), by = "sample_id") |>
  inner_join(
    select(chosen, uniprot_id, gene, selected_in, selection_p, facet),
    by = "uniprot_id"
  ) |>
  filter(!is.na(log2_abundance)) |>
  mutate(
    facet = factor(facet, chosen$facet),
    observed = corrected$other$n.observations[cbind(uniprot_id, sample_id)] > 0,
    value_source = factor(
      if_else(observed, "observed", "modelled (no precursor seen)"),
      c("observed", "modelled (no precursor seen)")
    )
  ) |>
  select(
    uniprot_id, gene, selected_in, selection_p, sample_id, subject, group, timepoint,
    log2_abundance, observed, facet, value_source
  )

stopifnot(
  n_distinct(panel_data$uniprot_id) == n_shown, all(chosen$pi_score < alpha_nominal),
  all(panel_data$log2_abundance == corrected$E[cbind(panel_data$uniprot_id, panel_data$sample_id)])
)

means <- summarise(panel_data,
  log2_abundance = mean(log2_abundance),
  .by = c(facet, group, timepoint)
)
panel <- panel_data |>
  ggplot(aes(timepoint, log2_abundance, colour = group)) +
  geom_line(aes(group = subject), alpha = 0.6, linewidth = 0.25) +
  geom_line(aes(group = group), means, linewidth = 0.9) +
  geom_point(aes(fill = value_source), size = 0.7, shape = 21, stroke = 0.25) +
  geom_point(aes(shape = group), means, size = 1.4) +
  facet_wrap(~facet, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = group_colours, name = NULL) +
  scale_shape_manual(values = group_shapes, name = NULL) +
  scale_fill_manual(
    values = c(observed = "grey40", `modelled (no precursor seen)` = "white"), name = NULL
  ) +
  scale_x_discrete(expand = expansion(add = 0.2)) +
  labs(
    title = "HR - LR candidates, person by person",
    subtitle = sprintf(
      paste0(
        "The %d smallest p among panel b's candidates;\nblood-corrected values. ",
        "Thin lines: people;\nthick: group means. Open points: values limpa\n",
        "modelled (no precursor seen); (n/%d): biopsies observed"
      ),
      n_shown, ncol(corrected$E)
    ),
    x = NULL, y = "log2 abundance"
  ) +
  figure_theme +
  theme(
    strip.text = element_text(size = 5.5, hjust = 0), legend.position = "bottom",
    legend.box = "vertical", legend.margin = margin(t = -4)
  )
panel_data <- select(panel_data, -facet, -value_source)

save_panel(panel, "F2", "c_trajectories", 75, 135)
