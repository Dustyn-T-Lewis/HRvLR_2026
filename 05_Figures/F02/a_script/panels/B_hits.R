# F02 B: nominal proteins per contrast against the 5% chance expects; none reaches BH < 0.05.
source(here::here("05_Figures", "theme.R"))

data <- read_sheet(
  "02_Differential_Expression/02_Differential/c_data/02_differential.xlsx", "contrast_summary"
) |>
  mutate(n_expected = 0.05 * n_tested, contrast = factor(contrast, rev(contrast_order)))

plot <- ggplot(data, aes(y = contrast)) +
  geom_col(aes(x = n_nominal), fill = "grey70", width = 0.7) +
  geom_point(aes(x = n_expected), shape = 124, size = 3) +
  labs(
    title = "Nominal proteins per contrast", x = "proteins at p < 0.05 (tick: 5% of tested)",
    y = NULL
  ) +
  theme_figure()
save_panel(plot, "F02", "B_hits")
