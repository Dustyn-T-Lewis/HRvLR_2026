# The protein dot-matrix (F2): one row per protein, labelled by gene, and one column per contrast
# under its window strip. Colour is logFC (PiYG), size -log10 p; a filled dot is a Π hit
# (p^|logFC| < alpha_nominal), an open one is not, and a black outline marks FDR < alpha_fdr.
# A filled dot carries logFC as its fill, an open one as its outline. Colour is capped (see below);
# the logFC in the panel data is not.
# Each row reads "gene (observed/n_samples)", the biopsies where a precursor was seen. A protein
# seen in few biopsies takes most of its values from limpa's detection model, which can give it a
# large logFC on little data.
# `results` has one row per chosen protein and contrast, with the protein_by_contrast columns plus
# `block` (row facet, a factor) and `order` (row position within the block, 1 at the top).
# `contrasts` gives the column order. Expects dplyr and ggplot2 attached and figure_style.R sourced.
pi_calls <- c(sprintf("Π < %s", alpha_nominal), sprintf("Π ≥ %s", alpha_nominal))
fdr_call <- sprintf("FDR < %s", alpha_fdr)

protein_grid <- function(results, contrasts, title, subtitle, n_samples, p_cap = 4) {
  d <- results |>
    mutate(
      family = contrast_panel(contrast),
      column = factor(contrast, contrasts),
      row = reorder(paste(block, uniprot_id), -order),
      log_p = pmin(-log10(p), p_cap),
      pi_call = factor(if_else(pi_score < alpha_nominal, pi_calls[1], pi_calls[2]), pi_calls)
    )
  # A gene shared by two accessions gets the accession, so no two rows read alike.
  labels <- distinct(d, row, gene, uniprot_id, n_observed) |>
    mutate(
      accession = if_else(
        duplicated(gene) | duplicated(gene, fromLast = TRUE), paste0(uniprot_id, ", "), ""
      ),
      label = sprintf("%s (%s%d/%d)", gene, accession, as.integer(n_observed), n_samples)
    )
  # A few proteins (blood neutrophil markers after the bout) move far more than the rest; the
  # colour stops at the median |logFC| of the drawn Π hits so the rest stay readable.
  limit <- ceiling(quantile(abs(d$logFC[d$pi_score < alpha_nominal]), 0.5) * 4) / 4
  strips <- contrast_strips(x = levels(d$family), y = levels(d$block))
  ggplot(d, aes(column, row)) +
    geom_point(
      aes(size = log_p, colour = logFC, shape = pi_call), filter(d, pi_call == pi_calls[2]),
      stroke = 0.4
    ) +
    geom_point(
      aes(size = log_p, fill = logFC, shape = pi_call), filter(d, pi_call == pi_calls[1]),
      colour = "grey30", stroke = 0.25
    ) +
    geom_point(
      aes(size = log_p, alpha = fdr_call), filter(d, fdr < alpha_fdr),
      shape = 1, colour = "black", stroke = 0.7
    ) +
    scale_colour_distiller(
      palette = "PiYG", direction = 1, limits = c(-limit, limit), oob = scales::squish,
      aesthetics = c("colour", "fill"),
      name = sprintf("logFC\n(capped at\n\u00b1%s)", limit)
    ) +
    scale_size_area(
      max_size = 2.6, limits = c(0, p_cap), breaks = c(1, 2, p_cap),
      labels = c("1", "2", paste0(p_cap, "+")), name = "-log10 p"
    ) +
    scale_shape_manual(values = rlang::set_names(c(21, 1), pi_calls), name = NULL, drop = FALSE) +
    scale_alpha_manual(values = rlang::set_names(1, fdr_call), name = NULL) +
    scale_x_discrete(labels = rlang::set_names(contrast_markup(contrasts), contrasts)) +
    scale_y_discrete(labels = tibble::deframe(select(labels, row, label))) +
    guides(
      colour = guide_colourbar(order = 1, barheight = unit(14, "mm"), barwidth = unit(2.5, "mm")),
      fill = "none",
      size = guide_legend(
        order = 2,
        override.aes = list(shape = 21, fill = "grey60", colour = "grey30", stroke = 0.25)
      ), shape = guide_legend(
        order = 3, override.aes = list(size = 2, fill = "grey60", colour = "grey30")
      ),
      alpha = guide_legend(order = 4, override.aes = list(size = 2))
    ) +
    ggh4x::facet_grid2(
      block ~ family,
      scales = "free", space = "free", switch = "y", strip = strips,
      labeller = labeller(family = label_wrap_gen(10))
    ) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    figure_theme +
    theme(
      axis.text.y = element_text(size = 5.2),
      axis.text.x.bottom = ggtext::element_markdown(size = 5.5, angle = 45, hjust = 1, vjust = 1),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90, size = 5.5),
      strip.text.x = element_text(size = 5.5),
      panel.spacing = unit(1, "mm"),
      legend.spacing.y = unit(1, "mm"),
      legend.text = element_text(size = 5.5), legend.title = element_text(size = 6)
    )
}
