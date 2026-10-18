# F1 c: each group's T1 to T2 change in fibre area, muscle area and strength, as a geometric-mean %
# with its 95% CI on a log-ratio axis, and the HR - LR p beside each trait. Muscle area and mixed
# fibre area formed the groups, so they sit apart and carry no HR - LR p.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(slopes = "02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.xlsx")
holds <- paste(
  "Panel c: T1 to T2 % change per trait and group with 95% CI (t-test on log ratios), and the",
  "HR - LR difference (Welch)."
)
effects <- read_sheet(inputs[["slopes"]], "trait_effects")

blocks <- c("formed the groups", "not used to form them")
difference <- filter(effects, group == "HR - LR")
panel_data <- effects |>
  filter(group != "HR - LR") |>
  left_join(select(difference, trait, hr_lr_pct = pct_change, hr_lr_p = p), by = "trait") |>
  mutate(
    block = factor(if_else(used_for_groups, blocks[1], blocks[2]), blocks),
    label = if_else(
      used_for_groups, trait_label(trait),
      sprintf("%s\nHR - LR p %s", trait_label(trait), signif(hr_lr_p, 2))
    ),
    trait = factor(trait, rev(trait_order)),
    group = factor(group, c("LR", "HR")),
    called = p < alpha_nominal
  ) |>
  arrange(block, desc(trait), desc(group))
stopifnot(
  nrow(panel_data) == 2 * nrow(difference), setequal(panel_data$group, c("HR", "LR")),
  all(panel_data$lower <= panel_data$pct_change & panel_data$pct_change <= panel_data$upper)
)

ratio <- \(pct) 1 + pct / 100
breaks <- c(0.75, 1, 1.5)
panel <- ggplot(panel_data, aes(ratio(pct_change), trait, colour = group)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey50") +
  geom_linerange(
    aes(xmin = ratio(lower), xmax = ratio(upper)),
    position = position_dodge(width = 0.6), linewidth = 0.45
  ) +
  geom_point(
    aes(shape = group, fill = if_else(called, group, NA)),
    position = position_dodge(width = 0.6), size = 1.6, stroke = 0.5
  ) +
  scale_x_continuous(
    transform = "log2", breaks = breaks, labels = \(x) sprintf("%+.0f%%", 100 * (x - 1))
  ) +
  scale_y_discrete(labels = tibble::deframe(distinct(panel_data, trait, label))) +
  scale_colour_manual(values = group_colours, breaks = c("HR", "LR"), name = NULL) +
  scale_fill_manual(values = group_colours, na.value = "white", guide = "none") +
  scale_shape_manual(values = c(HR = 21, LR = 24), breaks = c("HR", "LR"), name = NULL) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    x = "% change, T1 to T2", y = NULL, title = "Trait change per group",
    subtitle = sprintf(
      "geometric mean, 95%% CI, log scale\nfilled: change p < %s; open: p ≥ %s",
      alpha_nominal, alpha_nominal
    )
  ) +
  figure_theme +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 5.5, colour = "grey30"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 5.8, lineheight = 0.9),
    legend.position = "bottom"
  )

save_panel(panel, "F1", "c_trait_effects", 70, 95)
