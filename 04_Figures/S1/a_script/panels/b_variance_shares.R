library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(variance = "02_HRvLR/04_Variance/c_data/variance.xlsx")
shares <- read_sheet(inputs[["variance"]], "variance_shares")
holds <- "Panel b: per protein, the share of its variance from each component (variancePartition)."

components <- c(
  subject = "person", timepoint = "time point", group = "HR/LR", blood = "blood index",
  Residuals = "residual"
)
component_colours <- c(
  person = "grey45", `time point` = time_colours[["T2"]],
  `HR/LR` = contrast_palette$line[contrast_palette$family == "HR - LR"],
  `blood index` = "#B2182B", residual = "grey80"
)
panel_data <- shares |>
  tidyr::pivot_longer(-uniprot_id, names_to = "component", values_to = "share") |>
  mutate(component = factor(components[component], components))
medians <- summarise(panel_data, median = median(share), .by = component)
stopifnot(
  !anyNA(panel_data$component),
  all(abs(rowSums(select(shares, -uniprot_id)) - 1) < 1e-6)
)

panel <- ggplot(panel_data, aes(component, share)) +
  geom_violin(
    aes(fill = component),
    scale = "width", linewidth = 0.2, colour = "grey30", alpha = 0.6
  ) +
  geom_boxplot(width = 0.15, outliers = FALSE, linewidth = 0.25, fill = "white") +
  geom_text(
    aes(y = 1.06, label = sprintf("%.0f%%", 100 * median)), medians,
    size = 2, colour = "grey20"
  ) +
  scale_fill_manual(values = component_colours, guide = "none") +
  scale_y_continuous(
    labels = scales::label_percent(), breaks = seq(0, 1, 0.25), limits = c(0, 1.1)
  ) +
  labs(
    x = NULL, y = "share of the protein's variance", title = "Where the variance sits",
    subtitle = sprintf(
      "variancePartition, %s proteins; number: median share",
      scales::comma(nrow(shares))
    )
  ) +
  figure_theme +
  theme(panel.grid.major.x = element_blank())

save_panel(panel, "S1", "b_variance_shares", 70, 95)
