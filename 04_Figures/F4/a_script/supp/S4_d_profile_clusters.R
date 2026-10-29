# Proteins that share a pattern of trait links: the largest profile groups of
# 03_Continuous/03_Profile_Clusters, each drawn as its members' mean limma t in every window-trait
# screen, named by the group's top over-represented set where that set passes FDR < 0.05.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "module_labels.R"))

inputs <- c(profile_clusters = "03_Continuous/03_Profile_Clusters/c_data/profile_clusters.rds")
profiles <- readRDS(here(inputs[["profile_clusters"]]))
n_groups <- 12

shown <- profiles$profile_groups |>
  slice_max(n_proteins, n = n_groups, with_ties = FALSE) |>
  mutate(row = sprintf("%d (%d) %s", group, n_proteins, set_identity(top_set, top_set_fdr, 40)))
panel_data <- profiles$profile_means |>
  filter(group %in% shown$group) |>
  separate_wider_delim(screen, " ", names = c("window", "trait"), cols_remove = FALSE) |>
  left_join(select(shown, group, n_proteins, top_set, top_set_fdr), by = "group") |>
  select(group, n_proteins, top_set, top_set_fdr, screen, window, trait, mean_t)

member_counts <- count(profiles$profile_members, group, name = "n_members")
stopifnot(
  n_distinct(panel_data$group) == n_groups,
  all(count(panel_data, group)$n == n_distinct(profiles$profile_means$screen)),
  setequal(panel_data$window, names(screen_windows)),
  all(shown$n_proteins == member_counts$n_members[match(shown$group, member_counts$group)]),
  isTRUE(all.equal(
    panel_data$mean_t,
    profiles$profile_means$mean_t[match(
      paste(panel_data$group, panel_data$screen),
      paste(profiles$profile_means$group, profiles$profile_means$screen)
    )]
  ))
)

drawn <- panel_data |>
  mutate(
    column = factor(trait_label(trait, change = FALSE), trait_label(trait_order)),
    window = factor(screen_windows[window], screen_windows),
    row = factor(shown$row[match(group, shown$group)], rev(shown$row))
  )
# Mean t is capped so the many small means stay visible beside a few large ones.
t_limit <- 1.5
panel <- ggplot(drawn, aes(column, row, fill = mean_t)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  ggh4x::facet_grid2(
    . ~ window,
    scales = "free_x", space = "free_x", strip = contrast_strips(x = screen_family)
  ) +
  scale_fill_distiller(
    palette = "RdBu", limits = c(-t_limit, t_limit), oob = scales::squish, name = "mean t"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Proteins that share a pattern of trait links",
    subtitle = sprintf(
      paste(
        "%d largest of %d profile groups (%s of %s proteins), mean limma t per screen;",
        "row: group (proteins) top over-represented set"
      ),
      n_groups, nrow(profiles$profile_groups), scales::comma(sum(shown$n_proteins)),
      scales::comma(nrow(profiles$profile_members))
    ),
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_blank(), panel.spacing.x = unit(1, "mm"),
    legend.key.height = unit(5, "mm"), legend.key.width = unit(2.5, "mm")
  )

save_panel(panel, "S4", "d_profile_clusters", 250, 70)
