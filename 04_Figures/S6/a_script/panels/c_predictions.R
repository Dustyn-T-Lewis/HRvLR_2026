# S6 c: for the lead pairing and trait, each person's held-out prediction against the observed
# trait, per feature set, HR and LR marked.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "elastic_net.R"))

inputs <- c(elastic_net = elastic_net_input)
holds <- "Panel c: per person and feature set, observed trait and held-out prediction for the lead."
panel_data <- in_lead(elastic_net$predictions) |>
  select(pairing, trait, features, subject, group, observed, predicted)
q2 <- in_lead(elastic_net$summary) |>
  transmute(features, text = sprintf("Q2 %.2f, p %.3f", q2, p_perm))
stopifnot(all(count(panel_data, features)$n >= 15))

d <- panel_data |>
  left_join(q2, by = "features") |>
  mutate(
    facet = factor(
      sprintf("%s\n%s", feature_set_names[features], text),
      sprintf("%s\n%s", feature_set_names[q2$features], q2$text)
    ),
    group = factor(group, c("HR", "LR"))
  )
panel <- ggplot(d, aes(observed, predicted, colour = group, shape = group)) +
  geom_abline(linetype = "22", colour = "grey60") +
  geom_point(size = 1.6) +
  scale_colour_manual(values = group_colours, name = NULL) +
  scale_shape_manual(values = group_shapes, name = NULL) +
  facet_wrap(~facet, nrow = 1) +
  labs(
    title = "Held-out predictions, person by person",
    subtitle = paste(
      "Each person predicted by models that never saw them (mean of 5 repeats);",
      "dashed: perfect"
    ),
    x = paste(
      trait_label(sub("^pct_", "", lead$trait), change = FALSE),
      if (startsWith(lead$trait, "pct_")) "(% change)" else ""
    ),
    y = "predicted"
  ) +
  figure_theme +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
save_panel(panel, "S6", "c_predictions", 125, 70)
