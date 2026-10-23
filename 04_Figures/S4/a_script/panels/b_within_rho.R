# The screen that holds the strongest slope difference (smallest p across every feature type,
# window and trait), feature by feature: Spearman rho with the trait within HR against within LR,
# with the slope-difference calls marked. Within-group rho rests on 8 people and is descriptive.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(group_trait_slopes = "02_HRvLR/07_Group_Trait_Slopes/c_data/group_trait_slopes.rds")
slopes <- readRDS(here(inputs[["group_trait_slopes"]]))
holds <- paste(
  "Panel b: Spearman rho within HR and within LR for every feature in the screen holding the",
  "smallest slope-difference p, with that feature's slope difference."
)
n_label <- 6

top <- slice_min(slopes$group_difference, p, n = 1, with_ties = FALSE)
in_screen <- \(d) {
  filter(d, feature_type == top$feature_type, window == top$window, trait == top$trait)
}
panel_data <- in_screen(slopes$within_group) |>
  select(feature, label, group, rho) |>
  pivot_wider(names_from = group, values_from = rho, names_prefix = "rho_") |>
  inner_join(
    select(in_screen(slopes$group_difference), feature, slope_difference, t, p, fdr),
    by = "feature"
  ) |>
  filter(!is.na(rho_HR), !is.na(rho_LR)) |>
  mutate(evidence = evidence_level(p, fdr))

stopifnot(
  top$feature %in% panel_data$feature,
  all(abs(c(panel_data$rho_HR, panel_data$rho_LR)) <= 1),
  nrow(panel_data) <= nrow(in_screen(slopes$group_difference)),
  !anyDuplicated(panel_data$feature)
)

labelled <- panel_data |>
  filter(evidence == evidence_levels[1]) |>
  slice_min(p, n = n_label) |>
  arrange(desc(rho_LR), rho_HR) |>
  mutate(
    short = if (top$feature_type == "sets") pathway_label(label, 30) else label,
    label_x = 0.05, label_y = seq(0.95, by = -0.09, length.out = n())
  )
n_hits <- sum(panel_data$evidence == evidence_levels[1])
n_crossed <- sum(
  panel_data$evidence == evidence_levels[1] & sign(panel_data$rho_HR) != sign(panel_data$rho_LR)
)
hit_colours <- setNames(
  c(contrast_tone("ink")[["Interaction"]], contrast_tone("ink")[["HR - LR"]], "grey80"),
  evidence_levels
)

panel <- ggplot(panel_data, aes(rho_HR, rho_LR)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_abline(linetype = "dashed", linewidth = 0.2, colour = "grey60") +
  geom_point(
    aes(colour = evidence),
    data = arrange(panel_data, desc(evidence)), size = 0.6
  ) +
  geom_segment(
    aes(xend = label_x, yend = label_y),
    labelled,
    linewidth = 0.2, colour = "grey40"
  ) +
  geom_label(
    aes(label_x, label_y, label = short), labelled,
    size = 1.8, hjust = 0, linewidth = 0, label.padding = unit(0.3, "mm"), fill = "white"
  ) +
  scale_colour_manual(values = hit_colours, name = "HR - LR slope", drop = FALSE) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    title = "The strongest screen, within each group",
    subtitle = sprintf(
      "%s, %s, %s %s; %d of %d FDR hits change sign",
      trait_label(top$trait), screen_windows[[top$window]], scales::comma(nrow(panel_data)),
      tolower(feature_types[[top$feature_type]]), n_crossed, n_hits
    ),
    x = "Spearman rho in HR", y = "Spearman rho in LR"
  ) +
  figure_theme +
  theme(legend.position = "bottom", legend.margin = margin(t = -4)) +
  guides(colour = guide_legend(override.aes = list(size = 1.5)))

ggsave(
  here("04_Figures", "S4", "b_reports", "panels", "b_within_rho.pdf"), panel,
  width = 85, height = 90, units = "mm", device = cairo_pdf
)
