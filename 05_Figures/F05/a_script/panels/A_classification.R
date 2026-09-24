# F05 A: nominal classification hits over chance per task, at each level.
source(here::here("05_Figures", "theme.R"))

tasks <- c(
  "Training, HR (T1 to T2)", "Training, LR (T1 to T2)", "Acute bout, HR (T2 to T3)",
  "Acute bout, LR (T2 to T3)", "HR vs LR at T1 (floor)", "HR vs LR at T2",
  "HR vs LR, training change", "HR vs LR, acute change"
)
data <- read_chance() |>
  filter(analysis == "classification") |>
  mutate(comparison = factor(comparison, rev(tasks)))

plot <- ggplot(data, aes(ratio, comparison, colour = level)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = level_colours, name = NULL) +
  labs(title = "Classification", x = "nominal hits / chance expectation", y = NULL) +
  theme_figure()
save_panel(plot, "F05", "A_classification")
