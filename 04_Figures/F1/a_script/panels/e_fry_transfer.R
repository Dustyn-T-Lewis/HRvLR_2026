# F1 e: each group's training and acute responders tested with fry in the other group's contrast.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx")
holds <- paste(
  "Panel e: each group's training and acute sets (Π or p < 0.05) tested with fry in the other",
  "group's contrast; sets under 5 proteins were not tested."
)
fry_rows <- read_sheet(inputs[["contrasts"]], "fry_concordance")

criteria <- c(pi = sprintf("Π < %s", alpha_nominal), nominal = sprintf("p < %s", alpha_nominal))
panel_data <- fry_rows |>
  mutate(
    window = sub("_.*", "", set_from), from_group = sub(".*_", "", set_from),
    in_group = sub(".*_", "", ranked_on),
    set_direction = if_else(set_direction == "up", "Up", "Down"),
    evidence = evidence_level(p, fdr),
    signed_log_p = if_else(set_direction == "Up", 1, -1) * -log10(p)
  ) |>
  select(
    window, from_group, in_group, criterion, set_direction, set_size, tested, set_from,
    ranked_on, direction, p, fdr, evidence, mean_logfc_own, mean_logfc_other, same_sign_share,
    signed_log_p
  ) |>
  arrange(window, from_group, criterion, set_direction)
stopifnot(
  nrow(panel_data) == nrow(fry_rows), all(panel_data$from_group != panel_data$in_group),
  all(panel_data$tested == !is.na(panel_data$p))
)

windows <- c("Training", "Acute")
row_label <- function(window, from, to, criterion) {
  sprintf(
    "%s in %s, %s",
    in_colour(from, group_shade(window, from)), in_colour(to, group_shade(window, to)),
    criteria[criterion]
  )
}
plotted <- panel_data |>
  mutate(
    window = factor(window, windows),
    row = factor(paste(window, from_group, criterion)),
    side = if_else(set_direction == "Up", 1, -1),
    row_text = row_label(as.character(window), from_group, in_group, criterion),
    fry_direction = coalesce(direction, set_direction),
    text = case_when(
      !tested ~ sprintf("%d, too few", set_size),
      direction != set_direction ~ sprintf("%d · %.2g, opposite", set_size, fdr),
      .default = sprintf("%d · %.2g", set_size, fdr)
    ),
    text_x = if_else(tested, signed_log_p + side * 0.12, side * 0.12),
    hjust = if_else(side > 0, 0, 1)
  )
row_order <- plotted |>
  distinct(window, from_group, criterion, row, row_text) |>
  arrange(window, from_group, desc(criterion))
plotted <- mutate(plotted, row = factor(row, rev(row_order$row)))
# Room beyond the longest bar on each side for its label.
limits <- c(
  min(plotted$signed_log_p, 0, na.rm = TRUE) - 1.8,
  max(plotted$signed_log_p, 0, na.rm = TRUE) + 2.4
)

panel <- ggplot(filter(plotted, tested), aes(signed_log_p, row)) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_vline(
    xintercept = c(-1, 1) * -log10(alpha_nominal), linewidth = 0.3, linetype = "22",
    colour = "grey55"
  ) +
  geom_col(
    aes(fill = fry_direction, alpha = evidence, colour = fry_direction),
    width = 0.65, linewidth = 0.25
  ) +
  geom_text(
    aes(text_x, row, label = text, hjust = hjust), plotted,
    size = 1.9, colour = "grey15"
  ) +
  scale_fill_manual(
    values = direction_colours, aesthetics = c("fill", "colour"), name = NULL
  ) +
  scale_alpha_manual(values = c(1, 0.45, 0.12), drop = FALSE, name = NULL) +
  scale_y_discrete(labels = tibble::deframe(select(row_order, row, row_text))) +
  scale_x_continuous(
    labels = abs, limits = limits, expand = expansion(0)
  ) +
  ggh4x::facet_grid2(
    window ~ .,
    scales = "free_y", space = "free_y", switch = "y",
    strip = contrast_strips(y = windows)
  ) +
  labs(
    x = "-log10 p in the other group (down | up sets)", y = NULL,
    title = "Does one group's response carry over?",
    subtitle = paste0(
      "fry: proteins one group changes tested in the other group's contrast\n",
      "label: set size · FDR; dashed: p = 0.05"
    )
  ) +
  figure_theme +
  theme(
    axis.text.y.left = ggtext::element_markdown(size = 6),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 6, face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom", legend.box = "vertical",
    legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0),
    panel.spacing.y = unit(1.5, "mm")
  ) +
  guides(
    fill = guide_legend(order = 1), colour = guide_legend(order = 1),
    alpha = guide_legend(order = 2, override.aes = list(fill = "grey30", colour = "grey30"))
  )

save_panel(panel, "F1", "e_fry_transfer", 105, 95)
