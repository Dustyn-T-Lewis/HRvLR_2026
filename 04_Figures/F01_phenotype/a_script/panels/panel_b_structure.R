# Panel B: how the ten phenotypes relate to each other. This is what decides
# how many independent questions the association sweep is really asking - three
# fibre-area measures correlating above 0.9 are close to one variable tested
# three times, and BH inside a cell cannot see that.
if (!exists("pheno")) {
  source(here::here("04_Figures", "F01_phenotype", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, tidyr, tibble, dplyr)

build_structure <- function(pheno, tag = "B") {
  m <- stats::cor(
    pheno[, PHENOTYPES],
    use = "pairwise.complete.obs", method = "spearman"
  )
  ord <- stats::hclust(stats::as.dist(1 - m))$order
  lab <- TRAIT_LABELS[PHENOTYPES][ord]

  d <- as_tibble(m, rownames = "a") |>
    pivot_longer(-"a", names_to = "b", values_to = "rho") |>
    mutate(
      a = factor(TRAIT_LABELS[.data$a], levels = lab),
      b = factor(TRAIT_LABELS[.data$b], levels = rev(lab))
    )

  p <- ggplot(d, aes(.data$a, .data$b, fill = .data$rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", .data$rho)),
      size = 1.9,
      colour = "grey20"
    ) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = c(-1, 1), name = "rho"
    ) +
    labs(
      title = "The phenotypes are not ten independent questions", tag = tag,
      subtitle = "Spearman correlation, clustered",
      x = NULL, y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank(),
      legend.key.width = unit(6, "pt")
    )

  list(plot = p, audit = as_tibble(m, rownames = "trait"))
}

b_structure <- build_structure(pheno)
save_panel(
  b_structure$plot, file.path(F01_RPT, "panels", "panel_b_structure"), 135, 120
)
F01_PANELS[["structure"]] <- b_structure$plot
F01_AUDIT[["phenotype_correlations"]] <- b_structure$audit
