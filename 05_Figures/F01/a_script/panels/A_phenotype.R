# F01 A: four of the ten adaptation measures, one point per subject, by arm.
source(here::here("05_Figures", "theme.R"))

labels <- c(
  comp_hypertrophy = "composite hypertrophy", d_fcsa_mixed = "change in mixed fCSA (µm²)",
  d_mcsa = "change in mCSA (cm²)", d_1rm_legpress = "change in leg-press 1RM (kg)"
)
data <- readr::read_csv(here("00_Input", "phenotype.csv"), show_col_types = FALSE) |>
  select(subject, arm, all_of(names(labels))) |>
  pivot_longer(-c(subject, arm), names_to = "outcome") |>
  mutate(outcome = factor(labels[outcome], labels))

plot <- ggplot(data, aes(arm, value, colour = arm)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "grey40") +
  geom_point(position = position_jitter(width = 0.12, seed = 1), size = 1.2) +
  facet_wrap(~outcome, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = arm_colours, guide = "none") +
  labs(title = "Adaptation by arm", x = NULL, y = NULL) +
  theme_figure()
save_panel(plot, "F01", "A_phenotype", width = 178, height = 55)
