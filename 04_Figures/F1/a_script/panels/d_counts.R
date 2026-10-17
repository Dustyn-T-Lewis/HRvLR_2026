# F1 d: proteins per contrast at Π, p and FDR < 0.05, up and down, against the count chance gives.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx")
holds <- "Panel d: per contrast and direction, proteins at Π < 0.05, p < 0.05 and FDR < 0.05."
results <- read_sheet(inputs[["contrasts"]], "protein_by_contrast")
summary_rows <- read_sheet(inputs[["contrasts"]], "contrast_summary")

# Rows in the design's order, under their family's box: the group for a within-group window, the
# time point for HR - LR, the window for an interaction.
contrast_rows <- tibble(contrast = summary_rows$contrast) |>
  mutate(
    family = factor(contrast_family(contrast), c("Training", "Acute", "HR - LR", "Interaction")),
    label = contrast_short(contrast)
  ) |>
  arrange(family) |>
  mutate(contrast = factor(contrast, rev(contrast)))

panel_data <- results |>
  add_count(contrast, name = "n_tested") |>
  mutate(direction = factor(if_else(logFC > 0, "Up", "Down"), c("Up", "Down"))) |>
  summarise(
    n_pi = sum(sig_pi != 0), n_nominal = sum(p < alpha_nominal), n_fdr_05 = sum(fdr < alpha_fdr),
    n_tested = first(n_tested), .by = c(contrast, direction)
  ) |>
  complete(contrast, direction, fill = list(n_pi = 0, n_nominal = 0, n_fdr_05 = 0)) |>
  mutate(n_tested = max(n_tested, na.rm = TRUE), .by = contrast) |>
  inner_join(mutate(contrast_rows, contrast = as.character(contrast)), by = "contrast") |>
  mutate(n_chance = alpha_nominal * n_tested / 2) |>
  select(family, label, contrast, direction, n_tested, n_pi, n_nominal, n_fdr_05, n_chance) |>
  arrange(family, contrast, direction)

totals <- panel_data |>
  summarise(
    pi_up = sum(n_pi[direction == "Up"]), pi_down = sum(n_pi[direction == "Down"]),
    n_nominal = sum(n_nominal), n_fdr_05 = sum(n_fdr_05), n_tested = first(n_tested),
    .by = contrast
  ) |>
  inner_join(summary_rows, by = "contrast", suffix = c("", "_source"))
stopifnot(
  nrow(totals) == nrow(summary_rows),
  with(totals, all(
    pi_up == pi_up_source, pi_down == pi_down_source, n_nominal == n_nominal_source,
    n_fdr_05 == n_fdr_05_source, n_tested == n_tested_source
  ))
)

pi_text <- sprintf("Π < %s", alpha_nominal)
p_text <- sprintf("p < %s", alpha_nominal)
bar_fills <- rlang::set_names(
  c(direction_colours[["Down"]], "#F1B6DA", "#B8E186", direction_colours[["Up"]]),
  paste(c("Down", "Down", "Up", "Up"), c(pi_text, p_text, p_text, pi_text))
)
# Bar widths (the Π bar sits inside the p bar), the gap before an FDR count, and room on the left
# for the row boxes.
p_width <- 0.75
pi_width <- 0.38
fdr_gap <- 4
box_room <- 1.4
side <- \(direction) if_else(direction == "Up", 1, -1)
bars <- panel_data |>
  pivot_longer(c(n_nominal, n_pi), names_to = "tier", values_to = "n") |>
  mutate(
    signed = side(direction) * n,
    contrast = factor(contrast, levels(contrast_rows$contrast)),
    fill = paste(direction, if_else(tier == "n_pi", pi_text, p_text)),
    width = if_else(tier == "n_pi", pi_width, p_width)
  ) |>
  # The Π bar is drawn last, over the p bar.
  arrange(tier == "n_pi")
fdr_marks <- panel_data |>
  filter(n_fdr_05 > 0) |>
  mutate(
    contrast = factor(contrast, levels(contrast_rows$contrast)),
    x = side(direction) * (pmax(n_nominal, n_pi) + fdr_gap),
    hjust = if_else(direction == "Up", 0, 1)
  )
bar_limit <- max(panel_data$n_nominal, panel_data$n_chance)
row_boxes <- mutate(contrast_rows, x = -bar_limit * 1.08)
chance <- panel_data$n_chance[1]

panel <- ggplot(bars, aes(signed, contrast)) +
  geom_col(aes(fill = fill, width = width), position = "identity") +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_vline(xintercept = c(-chance, chance), linewidth = 0.3, linetype = "22", colour = "grey35") +
  geom_text(
    aes(x, contrast, label = n_fdr_05, hjust = hjust), fdr_marks,
    size = 2.2, fontface = "bold", colour = "black"
  ) +
  # Families share the fill scale with the bars; the legend lists the bars only.
  geom_label(
    aes(x, contrast, label = label, fill = family, colour = family), row_boxes,
    hjust = 1, size = 2.1, fontface = "bold", text.colour = "black",
    linewidth = 0.35, label.r = unit(0.6, "mm"), label.padding = unit(0.7, "mm")
  ) +
  scale_fill_manual(
    values = c(bar_fills, contrast_tone("fill")), breaks = names(bar_fills), name = NULL
  ) +
  scale_colour_manual(values = contrast_tone("line"), guide = "none") +
  scale_x_continuous(labels = abs, limits = c(-bar_limit * box_room, bar_limit * 1.08)) +
  facet_grid(family ~ ., scales = "free_y", space = "free_y") +
  coord_cartesian(clip = "off") +
  labs(
    x = "proteins (down | up)", y = NULL, title = "Proteins per contrast",
    subtitle = sprintf(
      paste0(
        "%s tested; dashed: chance at %s\n(5%% of tested, half each way)\n",
        "numbers: proteins at FDR < %s"
      ),
      scales::comma(panel_data$n_tested[1]), p_text, alpha_fdr
    )
  ) +
  figure_theme +
  theme(
    strip.text = element_blank(), axis.text.y = element_blank(),
    legend.position = "bottom", panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(1.2, "mm")
  ) +
  guides(fill = guide_legend(nrow = 2))

save_panel(panel, "F1", "d_counts", 50, 95)
