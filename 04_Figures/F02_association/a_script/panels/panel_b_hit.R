# Panel B: the one cell that cleared BH, drawn as the raw scatter it came from.
# A regression at n = 15 can clear a threshold on the strength of two or three
# subjects, and the only way to see whether that happened is to look at the
# points rather than the p-value.
if (!exists("survivors")) {
  source(here::here("04_Figures", "F02_association", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, ggrepel, dplyr)

build_hit <- function(survivors, pheno, tag = "B") {
  top <- survivors |> slice_min(.data$bh, n = 1, with_ties = FALSE)
  feat <- feature_matrices()[[top$level]]
  delta <- subject_window(feat, top$window)
  subjects <- colnames(delta)
  d <- tibble::tibble(
    subject = subjects,
    change = delta[top$feature, ],
    outcome = phenotype_vector(pheno, top$phenotype)[subjects]
  ) |>
    filter(!is.na(.data$outcome))

  rho <- stats::cor(d$change, d$outcome, method = "spearman")

  p <- ggplot(d, aes(.data$outcome, .data$change)) +
    geom_smooth(
      method = "lm", formula = y ~ x, colour = "#2166AC", fill = "grey85",
      linewidth = 0.6
    ) +
    geom_point(size = 2.4, colour = "grey25") +
    ggrepel::geom_text_repel(aes(label = .data$subject),
      size = 1.9,
      max.overlaps = 20, seed = 42
    ) +
    labs(
      title = sprintf("%s vs %s", top$feature, TRAIT_LABELS[top$phenotype]),
      tag = tag,
      subtitle = sprintf(
        "%s, n = %d, BH = %.4f, Spearman rho = %.2f",
        WINDOW_LABELS[top$window], nrow(d), top$bh, rho
      ),
      x = TRAIT_LABELS[top$phenotype],
      y = sprintf("%s %s", top$feature, window_family(top$window))
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(
      size = 6.5, face = "italic", colour = "grey40"
    ))

  list(plot = p, audit = d |> mutate(
    level = top$level, window = top$window, phenotype = top$phenotype,
    feature = top$feature, spearman = rho, .before = 1
  ))
}

b_hit <- build_hit(survivors, pheno)
save_panel(b_hit$plot, file.path(F02_RPT, "panels", "panel_b_hit"), 110, 100)
F02_PANELS[["hit"]] <- b_hit$plot
F02_AUDIT[["top_hit_points"]] <- b_hit$audit
