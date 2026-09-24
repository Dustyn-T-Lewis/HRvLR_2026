# F04 D: moderated t of each eigengene on the nine contrasts; labels mark nominal p.
source(here::here("05_Figures", "theme.R"))

modules <- read_sheet("04_Network/01_build_modules/c_data/01_build_modules.xlsx", "module_summary")
data <- read_sheet("04_Network/04_test_modules/c_data/04_test_modules.xlsx", "eigengene_tests") |>
  mutate(
    module = factor(module, rev(modules$module)), contrast = factor(contrast, contrast_order)
  )

plot <- ggplot(data, aes(contrast, module, fill = t)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = if_else(p < 0.05, sprintf("%.3f", p), "")), size = 1.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "t") +
  labs(title = "Eigengenes on the nine contrasts", x = NULL, y = NULL) +
  theme_figure() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(plot, "F04", "D_eigengene_tests")
