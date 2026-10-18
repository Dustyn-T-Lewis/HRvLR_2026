# PCA panels of 02_HRvLR/01_Protein_Contrasts (F1b, S1e): one plot per window (a time point's
# levels or a person's change), people in HR and LR colours inside each group's 95% ellipse, each
# carrying the window's PERMANOVA for group, gap k and forced-split adjusted Rand index. Expects
# dplyr, ggplot2 and patchwork attached and figure_style.R sourced.

pca_inputs <- c(
  contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx",
  variance = "02_HRvLR/04_Variance/c_data/variance.xlsx"
)
pca_titles <- c(
  T1 = "T1 baseline", T2 = "T2 trained", T3 = "T3 after the bout",
  Training = "Δ training (T2 − T1)", Acute = "Δ acute (T3 − T2)"
)
pca_title_colours <- c(
  time_colours, Training = contrast_tone("line")[["Training"]],
  Acute = contrast_tone("line")[["Acute"]]
)

# The PCA scores of `windows` with each window's tests, checked against the sheets they come from.
pca_data <- function(windows) {
  pca <- read_sheet(pca_inputs[["contrasts"]], "pca")
  tests <- read_sheet(pca_inputs[["contrasts"]], "permanova") |>
    filter(term == "group") |>
    select(window, permanova_n = n_samples, permanova_r2 = r2, permanova_p = p) |>
    inner_join(
      select(read_sheet(pca_inputs[["variance"]], "separation"), window, k_gap, ari_group),
      by = "window"
    )
  d <- pca |>
    filter(window %in% windows) |>
    inner_join(tests, by = "window") |>
    select(window, subject, group, PC1, PC2, pc1_share, pc2_share, permanova_n:ari_group) |>
    arrange(match(window, windows), group, subject)
  stopifnot(
    setequal(d$window, windows), all(count(d, window)$n == d$permanova_n[!duplicated(d$window)])
  )
  d
}

pca_plot <- function(d) {
  window <- d$window[1]
  centres <- summarise(d, across(c(PC1, PC2), mean), .by = group)
  tests <- with(d[1, ], sprintf(
    "R² %.2f, %s\ngap k = %d; ARI %.2f", permanova_r2, label_p(permanova_p), k_gap, ari_group
  ))
  axis_title <- \(i, share) sprintf("PC%d (%.0f%%)", i, 100 * share)
  ggplot(d, aes(PC1, PC2, colour = group)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey88") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey88") +
    stat_ellipse(
      aes(fill = group), geom = "polygon", type = "norm", level = 0.95, alpha = 0.15,
      linewidth = 0.35
    ) +
    geom_point(aes(shape = group), size = 1.6) +
    geom_point(data = centres, shape = 3, size = 1.8, stroke = 0.6) +
    # The two centres nearly coincide, so each label is pushed clear with a line to its centre.
    ggrepel::geom_label_repel(
      aes(label = group), centres,
      size = 2.2, fontface = "bold", fill = alpha("white", 0.85), label.size = 0,
      label.padding = unit(0.4, "mm"), box.padding = unit(4, "mm"), min.segment.length = 0,
      segment.size = 0.25, seed = 1
    ) +
    scale_colour_manual(values = group_colours, aesthetics = c("colour", "fill"), guide = "none") +
    scale_shape_manual(values = group_shapes, guide = "none") +
    labs(
      title = pca_titles[[window]], subtitle = tests,
      x = axis_title(1, d$pc1_share[1]), y = axis_title(2, d$pc2_share[1])
    ) +
    figure_theme +
    theme(
      plot.title = element_text(colour = pca_title_colours[[window]], face = "bold", size = 6.5),
      plot.subtitle = element_text(size = 5.5, colour = "grey25", lineheight = 0.95),
      axis.title = element_text(size = 6)
    )
}

# One plot per window in a row, wrapped so a composite tags the row once.
pca_row <- function(d, title, subtitle) {
  plots <- lapply(split(d, factor(d$window, unique(d$window))), pca_plot)
  wrap_elements(panel = wrap_plots(plots, nrow = 1)) +
    labs(title = title, subtitle = subtitle) +
    figure_theme
}
