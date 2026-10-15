library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx")
read_sheet <- \(sheet) readxl::read_excel(here(inputs[["contrasts"]]), sheet)
fry_rows <- read_sheet("fry_concordance")
summary_rows <- read_sheet("contrast_summary")
results <- read_sheet("protein_by_contrast")
holds <- paste(
  "Panel d: each group's training and acute sets (Π or p < 0.05) tested with fry in the other",
  "group's contrast; sets under 5 proteins were not tested."
)

# Every set the step could have formed, so a set too small for fry still gets its row.
criteria <- c(pi = "Π < 0.05", nominal = "p < 0.05")
set_sizes <- results |>
  filter(contrast %in% c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR")) |>
  mutate(set_direction = if_else(logFC > 0, "Up", "Down")) |>
  summarise(
    pi = sum(sig_pi != 0), nominal = sum(p < alpha_nominal),
    .by = c(contrast, set_direction)
  ) |>
  pivot_longer(c(pi, nominal), names_to = "criterion", values_to = "set_size") |>
  rename(set_from = contrast)

panel_data <- set_sizes |>
  mutate(
    set = paste(criterion, tolower(set_direction), sep = "_"),
    window = sub("_.*", "", set_from), from_group = sub(".*_", "", set_from),
    in_group = if_else(from_group == "HR", "LR", "HR")
  ) |>
  left_join(select(fry_rows, -n_proteins), by = c("set_from", "set")) |>
  mutate(
    tested = !is.na(p),
    evidence = evidence_level(p, fdr),
    signed_log_p = if_else(set_direction == "Up", 1, -1) * -log10(p)
  ) |>
  select(
    window, from_group, in_group, criterion, set_direction, set_size, tested, set_from,
    ranked_on, direction, p, fdr, evidence, mean_logfc_own, mean_logfc_other, same_sign_share,
    signed_log_p
  ) |>
  arrange(window, from_group, criterion, set_direction)

traced <- inner_join(
  filter(panel_data, tested), fry_rows,
  by = c("set_from", "ranked_on", "direction", "p", "fdr")
)
stopifnot(
  nrow(traced) == nrow(fry_rows), all(traced$set_size == traced$n_proteins),
  all(panel_data$tested == (panel_data$set_size >= 5))
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
      direction != set_direction ~ sprintf("%d, FDR %.2g, opposite", set_size, fdr),
      .default = sprintf("%d, FDR %.2g", set_size, fdr)
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
  min(plotted$signed_log_p, 0, na.rm = TRUE) - 2.2, max(plotted$signed_log_p, 0, na.rm = TRUE) + 2.8
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
    x = "-log10 p in the other group (down sets left, up sets right)", y = NULL,
    title = "Does one group's response carry over to the other?",
    subtitle = paste0(
      "fry: proteins one group changes (Π or p < 0.05) tested in the other group's contrast\n",
      "label: set size and FDR; dashed: p = 0.05"
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

ggsave(
  here("04_Figures", "F1", "b_reports", "panels", "d_fry_transfer.pdf"), panel,
  width = 95, height = 85, units = "mm", device = cairo_pdf
)
