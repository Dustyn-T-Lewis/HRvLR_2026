# F02 C: p-value histograms for the primary contrast and the floor.
source(here::here("05_Figures", "theme.R"))

shown <- c("Training_Interaction", "Baseline_HRvLR")
data <- read_sheet(
  "02_Differential_Expression/02_Differential/c_data/02_differential.xlsx", "DEP_matrix"
) |>
  select(uniprot_id, paste0(shown, "_p")) |>
  pivot_longer(-uniprot_id, names_to = "contrast", values_to = "p") |>
  filter(!is.na(p)) |>
  mutate(contrast = factor(sub("_p$", "", contrast), shown))

plot <- ggplot(data, aes(p)) +
  geom_histogram(binwidth = 0.05, boundary = 0, fill = "grey45") +
  geom_hline(
    data = count(data, contrast) |> mutate(flat = n / 20), aes(yintercept = flat),
    linetype = 2, colour = "firebrick"
  ) +
  facet_wrap(~contrast, ncol = 1) +
  labs(title = "p-value distribution", x = "p", y = "proteins") +
  theme_figure()
save_panel(plot, "F02", "C_p_histogram")
