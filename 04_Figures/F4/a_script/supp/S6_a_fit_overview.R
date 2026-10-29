# S6 a: every elastic-net fit's held-out Q2, by pairing, trait and feature set, outlined where the
# permutation p is below 0.05 (solid) or 0.10 (dashed).
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "elastic_net.R"))

panel_data <- select(
  elastic_net$summary, pairing, trait, features, n_people, n_features, q2, null_95, p_perm,
  n_stable, n_recurring
)
stopifnot(nrow(panel_data) == 81, !anyDuplicated(select(panel_data, pairing, trait, features)))

d <- panel_data |>
  mutate(
    pairing = factor(pairing, unique(elastic_net$summary$pairing)),
    trait = factor(trait_label(sub("^pct_", "", trait), change = FALSE), unique(
      trait_label(sub("^pct_", "", elastic_net$summary$trait), change = FALSE)
    )),
    features = factor(feature_set_names[features], feature_set_names),
    lead = case_when(p_perm < 0.05 ~ "p < 0.05", p_perm < 0.10 ~ "p < 0.10", .default = NA)
  )
panel <- ggplot(d, aes(trait, pairing, fill = pmax(q2, -0.5))) +
  geom_tile(colour = "white") +
  geom_tile(
    aes(linetype = lead), filter(d, !is.na(lead)),
    fill = NA, colour = "black", linewidth = 0.5, width = 0.9, height = 0.9
  ) +
  geom_text(aes(label = sprintf("%.2f", q2)), size = 1.8) +
  scale_fill_gradient2(
    low = "grey85", mid = "white", high = "#1B7837", midpoint = 0, limits = c(-0.5, 0.4),
    name = "held-out Q2\n(floored at -0.5)"
  ) +
  scale_linetype_manual(
    values = c(`p < 0.05` = "solid", `p < 0.10` = "22"), name = "permutation",
    na.translate = FALSE
  ) +
  facet_grid(. ~ features) +
  labs(
    title = "Which feature sets carry a trait, everyone pooled",
    subtitle = sprintf(
      paste(
        "Elastic net, 16 people; Q2 from repeated 4-fold CV against 1,000 shuffles;",
        "%d of %d fits at p < 0.10 (about %d by chance)"
      ),
      sum(d$p_perm < 0.10), nrow(d), round(0.10 * nrow(d))
    ),
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank())
save_panel(panel, "S6", "a_fit_overview", 250, 80)
