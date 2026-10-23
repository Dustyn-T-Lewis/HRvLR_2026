library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(variance = "02_HRvLR/04_Variance/c_data/variance.xlsx")
pc_covariates <- read_sheet(inputs[["variance"]], "pc_covariates")
pc_scores <- read_sheet(inputs[["variance"]], "pc_scores")
holds <- "Panel d: R² of each covariate alone on each of PCs 1 to 5, and each PC's variance share."

covariates <- c(
  timepoint = "time point", group = "HR/LR", subject = "person", blood = "blood index",
  proteins_detected = "proteins detected", weight = "sample weight"
)
components <- unique(pc_covariates$component)
panel_data <- pc_covariates |>
  mutate(
    covariate_label = factor(covariates[covariate], rev(covariates)),
    component_label = sprintf("%s\n%.1f%%", component, 100 * pc_share),
    component_label = factor(component_label, unique(component_label))
  )
stopifnot(
  !anyNA(panel_data$covariate_label), nrow(panel_data) == length(components) * length(covariates),
  all(components %in% names(pc_scores))
)

panel <- ggplot(panel_data, aes(component_label, covariate_label, fill = r2)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.2f", r2), colour = r2 > 0.5),
    size = 2, show.legend = FALSE
  ) +
  scale_fill_distiller(
    palette = "Blues", direction = 1, limits = c(0, 1), name = "R²"
  ) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL, title = "What drives the main components",
    subtitle = "R² of each covariate alone; PC share of variance under its name"
  ) +
  figure_theme +
  theme(
    panel.grid.major = element_blank(), legend.position = "bottom",
    legend.title = element_text(vjust = 0.9),
    legend.key.width = unit(6, "mm"), legend.key.height = unit(2, "mm")
  )

save_panel(panel, "S1", "d_pc_covariates", 70, 70)
