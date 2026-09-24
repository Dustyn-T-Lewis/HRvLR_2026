# F01 B: proteins left after each filtering step.
source(here::here("05_Figures", "theme.R"))

data <- read_sheet("01_Preprocess/01_Filtering/c_data/01_filter.xlsx", "filter_log") |>
  mutate(step = sub("^remove: ", "", step), step = factor(step, rev(step)))

plot <- ggplot(data, aes(n_after, step)) +
  geom_col(fill = "grey70", width = 0.7) +
  geom_text(aes(label = n_after), hjust = -0.15, size = 2.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Proteins after each filter step", x = "proteins", y = NULL) +
  theme_figure()
save_panel(plot, "F01", "B_filtering")
