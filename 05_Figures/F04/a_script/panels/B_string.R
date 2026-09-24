# F04 B: STRING edges over expectation per module, labelled with its top set at FDR < 0.05.
source(here::here("05_Figures", "theme.R"))

data <- read_sheet(
  "04_Network/02_characterise_modules/c_data/02_characterise_modules.xlsx", "labels"
) |>
  mutate(
    top = if_else(top_set_fdr < 0.05, enrichVolcano::ev_clean_label(top_set), ""),
    top = gsub("\n", " ", top),
    module = factor(module, module[order(string_ratio)])
  )

plot <- ggplot(data, aes(string_ratio, module, fill = module)) +
  geom_col(colour = "grey30", linewidth = 0.2, width = 0.75) +
  geom_text(aes(label = top), hjust = -0.05, size = 2) +
  scale_fill_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.9))) +
  labs(title = "STRING edges and top set", x = "observed / expected edges", y = NULL) +
  theme_figure()
save_panel(plot, "F04", "B_string")
