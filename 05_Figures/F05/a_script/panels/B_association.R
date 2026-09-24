# F05 B: nominal phenotype associations over chance per window and outcome, at each level.
source(here::here("05_Figures", "theme.R"))

data <- read_chance() |>
  filter(analysis != "classification") |>
  mutate(window = factor(sub("association: ", "", analysis), c("training", "baseline", "acute")))

plot <- ggplot(data, aes(ratio, comparison, colour = level)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
  facet_wrap(~window, nrow = 1) +
  scale_colour_manual(values = level_colours, name = NULL) +
  labs(title = "Association with phenotype", x = "nominal hits / chance expectation", y = NULL) +
  theme_figure()
save_panel(plot, "F05", "B_association", width = 178)
