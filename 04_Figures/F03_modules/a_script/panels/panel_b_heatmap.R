# Panel B: every module against every phenotype, in every window. 720 tests on
# one grid. Marking only the cells that clear BH inside their own cell, because
# none of them clears the correction across the sweep and drawing them as
# findings would misstate the result.
if (!exists("module_results")) {
  source(here::here("04_Figures", "F03_modules", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats, dplyr)

build_heatmap <- function(module_results, tag = "B") {
  d <- module_results |>
    mutate(
      module = sub("^ME_", "", .data$feature),
      trait = factor(TRAIT_LABELS[.data$phenotype],
        levels = rev(unname(TRAIT_LABELS))
      ),
      win = factor(WINDOW_LABELS[.data$window],
        levels = unname(WINDOW_LABELS)
      ),
      label = fct_rev(factor(module_label(.data$module),
        levels = sort(unique(module_label(.data$module)))
      )),
      signed = -log10(.data$p) * sign(.data$slope),
      hit = .data$bh < 0.05
    )

  p <- ggplot(d, aes(.data$win, .data$label, fill = .data$signed)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_point(
      data = filter(d, .data$hit), shape = 4, size = 1.4, stroke = 0.6,
      colour = "black"
    ) +
    facet_wrap(~trait, nrow = 2) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
      name = expression(-log[10] * "(p) x sign")
    ) +
    labs(
      title = "Every module, phenotype and window", tag = tag,
      subtitle = paste(
        "x marks BH < 0.05 within a cell.",
        "None survives correction across the 180-cell sweep"
      ),
      x = NULL, y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      axis.text.x = element_text(size = 5.5, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 5),
      strip.text = element_text(size = 6),
      panel.grid = element_blank(),
      legend.key.width = unit(5, "pt"),
      legend.title = element_text(size = 6)
    )

  list(plot = p, audit = d |> dplyr::select(
    module, window, phenotype, n, slope, t, p, bh
  ))
}

b_heatmap <- build_heatmap(module_results)
save_panel(
  b_heatmap$plot, file.path(F03_RPT, "panels", "panel_b_heatmap"), 210, 130
)
F03_PANELS[["heatmap"]] <- b_heatmap$plot
F03_AUDIT[["module_phenotype_grid"]] <- b_heatmap$audit
