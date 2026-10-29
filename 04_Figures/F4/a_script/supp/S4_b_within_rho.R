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
  mutate(short = if (top$feature_type == "sets") pathway_label(label, 30) else label)
n_hits <- sum(panel_data$evidence == evidence_levels[1])
n_crossed <- sum(
  panel_data$evidence == evidence_levels[1] & sign(panel_data$rho_HR) != sign(panel_data$rho_LR)
)
hit_colours <- rlang::set_names(
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
  ggrepel::geom_label_repel(
    aes(label = short), labelled,
    size = 1.8, linewidth = 0, label.padding = unit(0.3, "mm"), fill = "white",
    segment.size = 0.2, segment.colour = "grey40", min.segment.length = 0, seed = 1,
    direction = "y", nudge_x = 1.2, xlim = c(0.05, 1), ylim = c(0, 1),
    force = 10, box.padding = 0.35, max.overlaps = Inf, max.iter = 1e5, max.time = 3
  ) +
  scale_colour_manual(values = hit_colours, name = "HR - LR slope", drop = FALSE) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
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

save_panel(panel, "S4", "b_within_rho", 125, 95)
