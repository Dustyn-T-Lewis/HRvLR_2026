# S6 d: does the lead's link hold in both groups? For each feature selected for the lead, Spearman
# rho with the trait within HR and within LR, and the HR - LR slope p
# (02_HRvLR/07_Group_Trait_Slopes).
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "elastic_net.R"))

inputs <- c(
  elastic_net = elastic_net_input,
  group_trait_slopes = "02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds"
)
# The step-07 window that pairs the same feature value with the same trait.
slope_window <- c(
  "T1 level vs T1 trait" = "baseline_level", "T1 level predicts gain" = "baseline",
  "training change tracks gain" = "training", "T2 level vs T2 trait" = "post_level",
  "acute change vs gain" = "acute"
)
window <- slope_window[[lead$pairing]]
slopes <- readRDS(here(inputs[["group_trait_slopes"]]))
type_of <- c(proteins = "proteins", modules = "modules", hallmark = "sets")
selected <- in_lead(elastic_net$selection) |>
  transmute(feature_type = type_of[features], features, feature, label, frequency)
panel_data <- selected |>
  inner_join(
    filter(slopes$within_group, window == .env$window, trait == lead$trait) |>
      select(feature_type, feature, group, rho),
    by = c("feature_type", "feature")
  ) |>
  tidyr::pivot_wider(names_from = group, values_from = rho, names_prefix = "rho_") |>
  left_join(
    filter(slopes$group_difference, window == .env$window, trait == lead$trait) |>
      select(feature_type, feature, slope_p = p),
    by = c("feature_type", "feature")
  )
stopifnot(nrow(panel_data) == nrow(selected), !anyNA(c(panel_data$rho_HR, panel_data$rho_LR)))

d <- panel_data |>
  mutate(
    features = factor(feature_set_names[features], feature_set_names),
    label = reorder(paste(features, label), frequency),
    slope_text = sprintf("%.2f", slope_p)
  )
panel <- ggplot(d, aes(y = label)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_segment(aes(x = rho_LR, xend = rho_HR, yend = label), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(x = rho_LR, colour = "LR", shape = "LR"), size = 1.6) +
  geom_point(aes(x = rho_HR, colour = "HR", shape = "HR"), size = 1.6) +
  geom_text(
    aes(x = 1.08, label = slope_text, fontface = if_else(slope_p < 0.05, "bold", "plain")),
    size = 1.8, hjust = 0
  ) +
  scale_colour_manual(values = group_colours, name = NULL) +
  scale_shape_manual(values = group_shapes, name = NULL) +
  scale_y_discrete(labels = \(x) sub("^[^ ]+ (sets )?", "", x)) +
  scale_x_continuous(limits = c(-1, 1.25), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  ggh4x::facet_grid2(features ~ ., scales = "free_y", space = "free_y") +
  coord_cartesian(clip = "off") +
  labs(
    title = "Is the link carried by one group?",
    subtitle = paste(
      "Spearman rho with the trait within each group of 8;",
      "right: HR - LR slope p (bold < 0.05)"
    ),
    x = "Spearman rho", y = NULL
  ) +
  figure_theme +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, face = "bold"))
save_panel(panel, "S6", "d_hr_lr", 125, 110)
