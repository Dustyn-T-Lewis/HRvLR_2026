# F02 A: samples on the first two principal components after normalisation.
source(here::here("05_Figures", "theme.R"))

book <- "01_Preprocess/02_Normalization/c_data/02_normalize.xlsx"
variance <- read_sheet(book, "components")$variance_explained
data <- read_sheet(book, "pca_scores")

plot <- ggplot(data, aes(PC1, PC2, colour = timepoint, shape = arm)) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = timepoint_colours, name = NULL) +
  scale_shape_manual(values = c(HR = 16, LR = 1), name = NULL) +
  labs(
    title = "Samples after normalisation",
    x = sprintf("PC1 (%.1f%%)", 100 * variance[1]), y = sprintf("PC2 (%.1f%%)", 100 * variance[2])
  ) +
  theme_figure()
save_panel(plot, "F02", "A_pca")
