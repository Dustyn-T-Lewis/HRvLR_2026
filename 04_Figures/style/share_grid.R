# The trait-screen tile grid (F4b, S4a): one tile per window and trait, one facet per feature type,
# filled by the share of features at p < 0.05 (white at the 5% chance gives) and printed with that
# share in percent, plus the count at FDR < 0.05 in bold below it where there is one. Expects
# dplyr and ggplot2 attached and figure_style.R sourced.
#
# `counts`: one row per feature_type, window and trait, with share_nominal and n_fdr_05.
feature_types <- c(proteins = "Proteins", sets = "Gene sets", modules = "Modules")

# Check the saved counts against a recount of the per-feature results (p and fdr per feature_type,
# window and trait), so the tiles cannot drift from the tests behind them.
check_counts <- function(counts, results) {
  recount <- results |>
    summarise(
      n_tested = sum(!is.na(p)), n_nominal = sum(p < alpha_nominal, na.rm = TRUE),
      n_fdr_05 = sum(fdr < alpha_fdr, na.rm = TRUE), .by = c(feature_type, window, trait)
    ) |>
    arrange(feature_type, window, trait)
  saved <- counts |>
    select(feature_type, window, trait, n_tested, n_nominal, n_fdr_05) |>
    arrange(feature_type, window, trait)
  stopifnot(identical(saved, recount))
}

share_grid <- function(counts, title, subtitle) {
  window_text <- rlang::set_names(
    if_else(
      is.na(screen_family), screen_windows,
      in_colour(screen_windows, contrast_tone("ink")[screen_family])
    ),
    names(screen_windows)
  )
  drawn <- counts |>
    mutate(
      row = factor(trait_label(trait, change = FALSE), rev(trait_label(trait_order))),
      column = factor(window_text[window], window_text),
      feature_type = factor(feature_types[feature_type], feature_types),
      label = if_else(
        n_fdr_05 > 0, sprintf("%.0f\n%d", 100 * share_nominal, n_fdr_05),
        sprintf("%.0f", 100 * share_nominal)
      )
    )
  stopifnot(!anyNA(drawn$row), !anyNA(drawn$column), !anyNA(drawn$feature_type))
  ggplot(drawn, aes(column, row, fill = share_nominal)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(
      aes(label = label, fontface = if_else(n_fdr_05 > 0, "bold", "plain")),
      size = 1.8, lineheight = 0.85
    ) +
    facet_grid(. ~ feature_type) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D", midpoint = alpha_nominal,
      limits = c(0, 0.2), oob = scales::squish, labels = scales::label_percent(),
      name = sprintf("%% at\np < %s", alpha_nominal)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    figure_theme +
    theme(
      axis.text.x.bottom = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5),
      strip.text = element_text(face = "bold"), panel.grid.major = element_blank(),
      legend.key.height = unit(5, "mm"), legend.key.width = unit(2.5, "mm")
    )
}
