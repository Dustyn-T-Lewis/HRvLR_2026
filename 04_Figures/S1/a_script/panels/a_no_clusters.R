library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(variance = "02_HRvLR/04_Variance/c_data/variance.xlsx")
read_sheet <- \(sheet) readxl::read_excel(here(inputs[["variance"]]), sheet)
correlation <- read_sheet("biopsy_correlation")
biopsy_order <- read_sheet("biopsy_order")
gap <- read_sheet("gap")
biopsies <- read_sheet("pc_scores") |>
  select(sample_id, subject, group, timepoint)
holds <- paste(
  "Panel a: biopsy correlations in the saved clustering order, each biopsy's group, time point",
  "and person, and the gap statistic per k with the k chosen."
)

# The step's rule: the smallest k whose gap is within one standard error of the next k's.
k_chosen <- cluster::maxSE(gap$gap, gap$se, method = "Tibs2001SEmax")
strips <- biopsy_order |>
  inner_join(biopsies, by = "sample_id") |>
  arrange(position)
cells <- correlation |>
  inner_join(select(strips, sample_a = sample_id, position_a = position), by = "sample_a") |>
  inner_join(select(strips, sample_b = sample_id, position_b = position), by = "sample_b")
panel_data <- bind_rows(
  mutate(cells, layer = "correlation"),
  mutate(strips, layer = "strip"),
  mutate(gap, layer = "gap", chosen = k == k_chosen)
) |>
  relocate(layer)
n_biopsies <- nrow(strips)
stopifnot(
  nrow(cells) == n_biopsies^2, nrow(strips) == nrow(biopsies),
  setequal(strips$position, seq_len(n_biopsies)), all(gap$k == seq_len(nrow(gap)))
)

# The diagonal (r = 1) is left out of the colour range so the off-diagonal structure shows.
off_diagonal <- cells$r[cells$sample_a != cells$sample_b]
r_limit <- ceiling(max(abs(off_diagonal)) * 10) / 10
strip_rows <- c(group = -1.5, timepoint = -3.7, subject = -5.9)
person_colours <- rlang::set_names(
  grDevices::hcl.colors(n_distinct(strips$subject), "Dark 3"), unique(strips$subject)
)
strip_tiles <- strips |>
  tidyr::pivot_longer(c(group, timepoint, subject), names_to = "strip", values_to = "level") |>
  mutate(
    y = strip_rows[strip],
    colour = c(group_colours, time_colours, person_colours)[level]
  )

heatmap <- ggplot(cells, aes(position_a, n_biopsies + 1 - position_b)) +
  geom_raster(aes(fill = r)) +
  geom_tile(aes(position, y), strip_tiles, fill = strip_tiles$colour, height = 1.9) +
  annotate(
    "text",
    x = 0, y = strip_rows, label = c("HR/LR", "time", "person"),
    hjust = 1, size = 2, colour = "grey20"
  ) +
  scale_fill_distiller(
    palette = "RdBu", limits = c(-r_limit, r_limit), oob = scales::squish,
    name = "Pearson r"
  ) +
  scale_x_continuous(expand = expansion(add = c(6, 0.5))) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  coord_fixed(clip = "off") +
  labs(x = "biopsies, clustering order", y = NULL) +
  figure_theme +
  theme(
    axis.text = element_blank(), panel.grid.major = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(6, "mm"), legend.key.height = unit(2, "mm"),
    legend.title = element_text(vjust = 0.9)
  )

gap_plot <- ggplot(mutate(gap, chosen = k == k_chosen), aes(k, gap)) +
  geom_errorbar(aes(ymin = gap - se, ymax = gap + se), width = 0.25, linewidth = 0.3) +
  geom_line(linewidth = 0.3, colour = "grey50") +
  geom_point(aes(colour = chosen), size = 1.6, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey60")) +
  scale_x_continuous(breaks = gap$k) +
  labs(
    x = "k (clusters)", y = "gap ± SE", title = "Gap statistic",
    subtitle = sprintf("k = %d chosen", k_chosen)
  ) +
  figure_theme

legend_key <- tibble(
  label = c(names(group_colours), names(time_colours)),
  colour = c(group_colours, time_colours),
  x = c(1, 1, 2, 2, 2), y = c(2, 1, 3, 2, 1)
)
key_plot <- ggplot(legend_key, aes(x, y)) +
  geom_tile(fill = legend_key$colour, width = 0.18, height = 0.7) +
  geom_text(aes(x + 0.14, label = label), hjust = 0, size = 2) +
  scale_x_continuous(limits = c(0.85, 2.6)) +
  scale_y_continuous(limits = c(0.4, 3.6)) +
  theme_void()

# Wrapped as one element with its own title, so the composite tags the panel once.
grid <- (heatmap | (gap_plot / key_plot + plot_layout(heights = c(3, 1)))) +
  plot_layout(widths = c(3, 1))
panel <- wrap_elements(panel = grid) +
  labs(
    title = "Biopsies form no natural clusters",
    subtitle = sprintf(
      paste(
        "Correlation of %d biopsies, blood-corrected, in average-linkage order;",
        "%d people\nGap statistic, PAM on 10 PCs: k = %d"
      ),
      n_biopsies, n_distinct(strips$subject), k_chosen
    )
  ) +
  figure_theme

ggsave(
  here("04_Figures", "S1", "b_reports", "panels", "a_no_clusters.pdf"), panel,
  width = 110, height = 95, units = "mm", device = cairo_pdf
)
