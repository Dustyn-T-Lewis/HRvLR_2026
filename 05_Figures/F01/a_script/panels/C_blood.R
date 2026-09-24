# F01 C: the per-sample blood index on the 45 analysed samples, by arm and timepoint.
source(here::here("05_Figures", "theme.R"))

book <- "01_Preprocess/01_Filtering/c_data/01_filter.xlsx"
data <- read_sheet(book, "blood_index") |>
  inner_join(
    read_sheet(book, "outlier_diagnostics") |>
      filter(!consensus_outlier) |>
      select(sample_id, subject, arm, timepoint),
    by = "sample_id"
  )

plot <- ggplot(data, aes(timepoint, blood_index, colour = arm)) +
  geom_line(aes(group = subject), alpha = 0.3, linewidth = 0.3) +
  stat_summary(aes(group = arm), fun = mean, geom = "line", linewidth = 0.9) +
  geom_point(size = 0.9, alpha = 0.7) +
  scale_colour_manual(values = arm_colours, name = NULL) +
  labs(title = "Blood index by arm", x = NULL, y = "mean log2 haemoglobin intensity") +
  theme_figure()
save_panel(plot, "F01", "C_blood")
