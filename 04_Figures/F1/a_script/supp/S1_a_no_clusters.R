library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(variance = "02_HRvLR/04_Variance/c_data/variance.xlsx")
variance_sheet <- \(sheet) read_sheet(inputs[["variance"]], sheet)
correlation <- variance_sheet("biopsy_correlation")
biopsy_order <- variance_sheet("biopsy_order")
gap <- variance_sheet("gap")
biopsies <- variance_sheet("pc_scores") |>
  select(sample_id, subject, group, timepoint)
separation <- filter(variance_sheet("separation"), window == "all")

# The k the step chose (one-standard-error rule over the gap statistic) for all biopsies.
k_chosen <- separation$k_gap
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
  setequal(strips$position, seq_len(n_biopsies)), all(gap$k == seq_len(nrow(gap))),
  nrow(separation) == 1
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
  mutate(y = strip_rows[strip])

# Group and time point get a key; the person strip only shows the three biopsies of a person
# sitting together or apart.
strip_key <- c(group_colours, time_colours)
heatmap <- ggplot(cells, aes(position_a, n_biopsies + 1 - position_b)) +
  geom_raster(aes(fill = r)) +
  scale_fill_distiller(
    palette = "RdBu", limits = c(-r_limit, r_limit), oob = scales::squish,
    name = "Pearson r", guide = guide_colourbar(order = 1)
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(aes(position, y, fill = level), strip_tiles, height = 1.9) +
  scale_fill_manual(
    values = c(strip_key, person_colours), breaks = names(strip_key), name = NULL,
    guide = guide_legend(order = 2, nrow = 1)
  ) +
  annotate(
    "text",
    x = 0, y = strip_rows, label = c("HR/LR", "time", "person"),
    hjust = 1, size = 2, colour = "grey20"
  ) +
  scale_x_continuous(expand = expansion(add = c(6, 0.5))) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  coord_fixed(clip = "off") +
  labs(x = "biopsies, clustering order", y = NULL) +
  figure_theme +
  theme(
    axis.text = element_blank(), panel.grid.major = element_blank(),
    legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"),
    legend.key.width = unit(4, "mm"), legend.key.height = unit(2, "mm"),
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

# Wrapped as one element with its own title, so the composite tags the panel once.
grid <- (heatmap | gap_plot) + plot_layout(widths = c(3, 1))
panel <- wrap_elements(panel = grid) +
  labs(
    title = "Biopsies form no natural clusters",
    subtitle = sprintf(
      paste(
        "Correlation of %d biopsies, blood-corrected, in average-linkage order;",
        "%d people\nGap statistic, PAM on %d PCs: k = %d"
      ),
      n_biopsies, n_distinct(strips$subject), separation$n_pcs, k_chosen
    )
  ) +
  figure_theme

save_panel(panel, "S1", "a_no_clusters", 110, 95)
