# F04 C: Zsummary of each arm's modules tested in the other arm.
source(here::here("05_Figures", "theme.R"))

data <- read_sheet(
  "04_Network/03_preserve_modules/c_data/03_preserve_modules.xlsx", "preservation"
) |>
  mutate(direction = paste(reference, "modules in", test))

plot <- ggplot(data, aes(n_proteins, z_summary)) +
  geom_hline(yintercept = c(2, 10), linetype = 2, colour = "grey55") +
  geom_point(aes(fill = module), shape = 21, size = 1.8, colour = "grey30") +
  facet_wrap(~direction) +
  scale_fill_identity() +
  scale_x_log10() +
  labs(title = "Preservation between arms", x = "module size (log scale)", y = "Zsummary") +
  theme_figure()
save_panel(plot, "F04", "C_preservation")
