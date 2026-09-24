# F04 A: proteins per module and the subject ICC of its eigengene.
source(here::here("05_Figures", "theme.R"))

data <- read_sheet("04_Network/01_build_modules/c_data/01_build_modules.xlsx", "module_summary") |>
  mutate(module = factor(module, rev(module)))

plot <- ggplot(data, aes(n_proteins, module, fill = module)) +
  geom_col(colour = "grey30", linewidth = 0.2, width = 0.75) +
  geom_text(aes(label = sprintf("ICC %.2f", icc)), hjust = -0.15, size = 2) +
  scale_fill_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "Module sizes", x = "proteins", y = NULL) +
  theme_figure()
save_panel(plot, "F04", "A_modules")
