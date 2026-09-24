# F02 D: the primary contrast, every tested protein.
source(here::here("05_Figures", "theme.R"))

data <- read_sheet(
  "02_Differential_Expression/02_Differential/c_data/02_differential.xlsx", "DEP_matrix"
) |>
  transmute(
    uniprot_id, gene,
    log2fc = Training_Interaction_logFC, p = Training_Interaction_p,
    fdr = Training_Interaction_fdr
  ) |>
  filter(!is.na(p)) |>
  mutate(nominal = p < 0.05)

plot <- ggplot(data, aes(log2fc, -log10(p), colour = nominal)) +
  geom_point(size = 0.6, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "grey20"), guide = "none") +
  labs(
    title = "Training_Interaction", subtitle = sprintf("lowest FDR %.2f", min(data$fdr)),
    x = "log2 fold change", y = "-log10 p"
  ) +
  theme_figure()
save_panel(plot, "F02", "D_volcano")
