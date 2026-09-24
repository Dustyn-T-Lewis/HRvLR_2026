# F03 B: each set's HR against LR training NES; sets significant in either arm coloured.
source(here::here("05_Figures", "theme.R"))

data <- read_sheet(
  "03_Pathway_Enrichment/03_enrich_scatter_fgsea/c_data/03_enrich_scatter_fgsea.xlsx",
  "nes_scatter"
) |>
  filter(pair == "training") |>
  mutate(significance = factor(significance, c("Both", "Training_HR", "Training_LR", "NS")))

plot <- ggplot(data, aes(nes_x, nes_y)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point(data = filter(data, significance == "NS"), colour = "grey85", size = 0.4) +
  geom_point(data = filter(data, significance != "NS"), aes(colour = significance), size = 0.8) +
  scale_colour_manual(
    values = c(Both = "#6A3D9A", Training_HR = "#2166AC", Training_LR = "#B2182B"),
    name = "FDR < 0.05 in"
  ) +
  coord_fixed() +
  labs(title = "Training NES, HR against LR", x = "NES, Training_HR", y = "NES, Training_LR") +
  theme_figure()
save_panel(plot, "F03", "B_nes_training")
