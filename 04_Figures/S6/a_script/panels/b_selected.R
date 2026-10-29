# S6 b: for the lead pairing and trait, every feature kept in at least 40% of the 200 subsamples,
# per feature set, coloured by the direction of its link with the trait.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "elastic_net.R"))

inputs <- c(elastic_net = elastic_net_input)
holds <- "Panel b: per feature selected for the lead, its selection frequency, tier and rho."
panel_data <- in_lead(elastic_net$selection) |>
  select(pairing, trait, features, feature, label, frequency, tier, rho)
stopifnot(nrow(panel_data) > 0, all(panel_data$frequency >= 0.4))

d <- panel_data |>
  mutate(
    features = factor(feature_set_names[features], feature_set_names),
    label = reorder(paste(features, label), frequency),
    direction = if_else(rho > 0, "rises with the trait", "falls with the trait")
  )
panel <- ggplot(d, aes(frequency, label, fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0.4, linetype = "22", linewidth = 0.3) +
  geom_vline(xintercept = 0.6, linewidth = 0.3) +
  scale_fill_manual(
    values = c(`rises with the trait` = "#1B7837", `falls with the trait` = "#762A83"),
    name = NULL
  ) +
  scale_y_discrete(labels = \(x) sub("^[^ ]+ (sets )?", "", x)) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  ggh4x::facet_grid2(features ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = sprintf(
      "What carries it: %s, %s", lead$pairing,
      trait_label(sub("^pct_", "", lead$trait), change = FALSE)
    ),
    subtitle = paste(
      "Share of 200 subsamples of 12 people that keep the feature;",
      "solid: 0.6 (stable), dashed: 0.4 (recurring)"
    ),
    x = "selection frequency", y = NULL
  ) +
  figure_theme +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, face = "bold"))
save_panel(panel, "S6", "b_selected", 125, 110)
